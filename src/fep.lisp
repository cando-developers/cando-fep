(in-package :fep)



(defun core-group-p (group)
  (loop for atom in group
        for label = (chem:matter-get-property-or-default atom :label nil)
        when (and label (find #\? (string label)))
          do (return-from core-group-p t))
  nil)

(defun group-name (group)
  (loop for atom in group
        for label = (chem:matter-get-property-or-default atom :label nil)
        when label
          do (return-from group-name label)))

(defgeneric parse-residue-groups (sketch))

(defmethod parse-residue-groups ((sketch string))
  (multiple-value-bind (groups map-names-numbers)
      (jsme:parse-groups sketch)
    (let (core
          (side-chains (make-hash-table)))
      (loop for group in groups
            if (core-group-p group)
              do (setf core group)
            else do (setf (gethash (group-name group) side-chains) group))
      (values core side-chains map-names-numbers))))

(defmethod parse-residue-groups ((sketch chem:chem-draw))
  (let* ((agg (chem:as-aggregate sketch))
         core-molecule
         (core-atom-count 0)
         side-chain-molecules)
    (let ((num-bonds 0)
          (sum-lengths 0.0))
      (chem:map-bonds nil (lambda (a1 a2 order)
                            (let* ((v1 (chem:get-position a1))
                                   (v2 (chem:get-position a2))
                                   (delta (geom:v- v1 v2))
                                   (len (geom:vlength delta)))
                              (setf sum-lengths (+ sum-lengths len)
                                    num-bonds (1+ num-bonds))))
                      agg)
      (let* ((average-length (/ sum-lengths num-bonds))
             (scale-factor (/ 10.5 average-length)))
        (format t "Average bond length: ~a~%" average-length)
        ;; Scale bonds to 1.5 Angstroms
        (chem:map-atoms nil (lambda (a)
                          (let ((pos (chem:get-position a)))
                            (chem:set-position a (geom:v* pos scale-factor))))
                        agg)))
    (chem:map-molecules nil (lambda (m)
                              (let ((num-atoms (chem:number-of-atoms m)))
                                (when (> num-atoms core-atom-count)
                                  (setf core-molecule m
                                        core-atom-count num-atoms))))
                        agg)
    (chem:map-molecules nil (lambda (m)
                              (unless (eq m core-molecule)
                                (push m side-chain-molecules)))
                        agg)
    ;; Handle the core-group
    (let ((core-group (chem:map-atoms 'list #'identity core-molecule))
          (map-names-numbers (make-hash-table)))
      (mapcar (lambda (a)
                (if (eq (chem:get-element a) :ne)
                    (let ((bonds (chem:bonds-as-list a)))
                      (unless (= (length bonds) 1)
                        (warn "There are the wrong number of bonds from ~a bonds: ~a" a bonds))
                      (let* ((bond (first bonds))
                             (other-atom (chem:get-other-atom bond a))
                             (attach-sym (intern (string-trim "?" (string (chem:get-name a))) :keyword)))
                        (chem:set-property other-atom :attach attach-sym))
                      (chem:set-property a :label (chem:get-name a)))
                    (let* ((id (chem:matter-get-property-or-default a :id nil)))
                      (when id
                        (setf (gethash (chem:get-name a) map-names-numbers) id)))))
              core-group)
      ;; Handle the side-chains
      (let ((side-chain-groups (make-hash-table)))
        (loop for side-chain-molecule in side-chain-molecules
              do (let* (side-chain-name
                        (side-chain-atoms (chem:map-atoms
                                           'list
                                           (lambda (a)
                                             (when (eq (chem:get-element a) :ne)
                                               (setf side-chain-name (chem:get-name a))
                                               (let* ((other-name (intern (subseq (string side-chain-name) 0 2) :keyword))
                                                      (bonds (chem:bonds-as-list a)))
                                                 (unless (= (length bonds) 1)
                                                   (warn "There are the wrong number of bonds ~a bonds: ~a" a bonds))
                                                 (let* ((bond (first bonds))
                                                        (other-atom (chem:get-other-atom bond a)))
                                                   (chem:set-name other-atom other-name)
                                                   (chem:set-property a :label (chem:get-name a)))))
                                             a)
                                           side-chain-molecule)))
                   (setf (gethash side-chain-name side-chain-groups) side-chain-atoms)))
        (values core-group side-chain-groups map-names-numbers)))))



(defun labeled-atom (atom)
  (first (chem:bonded-atoms-as-list atom)))

(defclass bond-vector ()
  ((start-vector :initarg :start-vector :accessor start-vector)
   (direction :initarg :direction :accessor direction)))

(defun build-core-residue (atoms)
  (let ((residue (chem:make-residue :core))
        (out-atom-vectors (make-hash-table)))
    (loop for atom in atoms
          for label = (chem:matter-get-property-or-default atom :label nil)
          if label
            do (let* ((labeled-atom (labeled-atom atom))
                      (short-label (intern (subseq (string label) 0 2) :keyword))
                      (start (chem:get-position atom))
                      (dir (geom:v- (chem:get-position labeled-atom) start))
                      (bond-vector (make-instance 'bond-vector :start-vector start :direction dir)))
                 (chem:set-property labeled-atom :attach short-label)
                 (setf (gethash short-label out-atom-vectors) bond-vector)
                 (chem:remove-bond-to atom labeled-atom))
          else do (chem:add-matter residue atom))
    (values residue out-atom-vectors)))

(defun short-name (name)
  (intern (subseq (string name) 0 2) :keyword))

(defun atom-with-property (residue property property-value)
  (chem:map-atoms nil
                  (lambda (a)
                    (let ((pv (chem:matter-get-property-or-default a property nil)))
                      (when (eq pv property-value)
                        (return-from atom-with-property a))))
                  residue)
  nil)

(defun build-side-chain-residue (atoms side-chain-name core-residue core-vector)
  (let* ((residue (chem:make-residue side-chain-name))
         (short-side-chain-name (short-name side-chain-name))
         in-atom-vector attach-atom
         (core-atom (atom-with-property core-residue :attach short-side-chain-name)))
    (loop for atom in atoms
          for label = (chem:matter-get-property-or-default atom :label nil)
          if label
            do (let* ((labeled-atom (labeled-atom atom))
                      (short-label (intern (subseq (string label) 0 2) :keyword))
                      (start-vec (chem:get-position atom))
                      (dir (geom:v- (chem:get-position labeled-atom) start-vec)))
                 (setf attach-atom labeled-atom)
                 (chem:set-name labeled-atom short-label)
                 (setf in-atom-vector (make-instance 'bond-vector :start-vector start-vec :direction dir))
                 (chem:remove-bond-to atom labeled-atom))
          else do (chem:add-matter residue atom))
    (let* ((cross (geom:vcross (geom:vnormalized (direction in-atom-vector)) (geom:vnormalized (direction core-vector))))
           (angle (* 1.0 (asin (geom:vlength cross))))
           (residue-to-origin (geom:make-m4-translate (geom:v* (start-vector in-atom-vector) -1.0)))
           (residue-to-core (geom:make-m4-translate (chem:get-position core-atom)))
           (rotate (geom:make-m4-rotate-z angle)))
      (progn
        (chem:apply-transform-to-atoms residue residue-to-origin)
        (chem:apply-transform-to-atoms residue rotate)
        (chem:apply-transform-to-atoms residue residue-to-core)))
    (values residue in-atom-vector)))




(defvar *core-group* nil)
(defvar *side-chain-groups* nil)
(defvar *map-names-numbers* nil)

(defun build-topologys (sketch &key verbose)
  (let ((side-chains (make-hash-table))
        (topology-to-residue (make-hash-table)))
    (multiple-value-bind (core-group side-chain-groups map-names-numbers)
        (parse-residue-groups sketch)
      (setf *core-group* core-group
            *side-chain-groups* side-chain-groups
            *map-names-numbers* map-names-numbers)
      (format t "In build-topologys core-group: ~a~%" core-group)
      (format t "side-chain-groups: ~a~%" side-chain-groups)
      (format t "map-names-numbers: ~a~%" map-names-numbers)
      (multiple-value-bind (core-residue core-vectors)
          (build-core-residue core-group)
        (let ((core-topology (cando:make-simple-topology-from-residue core-residue)))
          (setf (gethash core-topology topology-to-residue) core-residue)
          (maphash (lambda (name vec)
                     (let* ((plug-name (intern (format nil "+~a" (string name)) :keyword))
                            (plug-atom (atom-with-property core-residue :attach name))
                            (plug-atom-name (chem:get-name plug-atom))
                            #+(or)(_ (format t "Making out plug plug-name: ~a    name: ~a~%" plug-name plug-atom-name))
                            (out-plug (chem:make-out-plug plug-name nil nil plug-atom-name :single-bond)))
                       (chem:add-plug core-topology (chem:get-name out-plug) out-plug)))
                   core-vectors)
          (loop for side-chain-name being the hash-keys in side-chain-groups using (hash-value side-chain-group)
                with side-chain-residue
                with side-chain-vector
                do (multiple-value-setq (side-chain-residue side-chain-vector)
                     (build-side-chain-residue side-chain-group
                                               side-chain-name
                                               core-residue
                                               (gethash (short-name side-chain-name) core-vectors)))
                do (let* ((side-chain-topology (cando:make-simple-topology-from-residue side-chain-residue))
                          (short-name (intern (subseq (string (chem:get-name side-chain-residue)) 0 2) :keyword))
                          (plug-name (intern (format nil "-~a" (string short-name)) :keyword))
                          (in-plug (chem:make-in-plug plug-name nil short-name :single-bond)))
                     (setf (gethash side-chain-topology topology-to-residue) side-chain-residue)
                     (chem:add-plug side-chain-topology plug-name in-plug)
                     (push (cons (chem:get-name side-chain-residue) side-chain-topology) (gethash short-name side-chains))))
          (values core-topology side-chains topology-to-residue map-names-numbers))))))



                                   
#+(or)
(defmethod build-topologys-chemdraw (chemdraw)
  (let ((side-chains (make-hash-table))
        (topology-to-residue (make-hash-table))
        (agg (chem:as-aggregate chemdraw)))
    (break "agg = ~a" agg)
    (multiple-value-bind (core-group side-chain-groups map-names-numbers)
        (parse-chemdraw chemdraw)
      (multiple-value-bind (core-residue core-vectors)
          (build-core-residue core-group)
        (let ((core-topology (cando:make-simple-topology-from-residue core-residue)))
          (setf (gethash core-topology topology-to-residue) core-residue)
          (maphash (lambda (name vec)
                     (let* ((plug-name (intern (format nil "+~a" (string name)) :keyword))
                            (plug-atom (atom-with-property core-residue :attach name))
                            (plug-atom-name (chem:get-name plug-atom))
                            #+(or)(_ (format t "Making out plug plug-name: ~a    name: ~a~%" plug-name plug-atom-name))
                            (out-plug (chem:make-out-plug plug-name nil nil plug-atom-name :single-bond)))
                       (chem:add-plug core-topology (chem:get-name out-plug) out-plug)))
                   core-vectors)
          (loop for side-chain-name being the hash-keys in side-chain-groups using (hash-value side-chain-group)
                with side-chain-residue
                with side-chain-vector
                do (multiple-value-setq (side-chain-residue side-chain-vector)
                     (build-side-chain-residue side-chain-group
                                               side-chain-name
                                               core-residue
                                               (gethash (short-name side-chain-name) core-vectors)))
                do (let* ((side-chain-topology (cando:make-simple-topology-from-residue side-chain-residue))
                          (short-name (intern (subseq (string (chem:get-name side-chain-residue)) 0 2) :keyword))
                          (plug-name (intern (format nil "-~a" (string short-name)) :keyword))
                          (in-plug (chem:make-in-plug plug-name nil short-name :single-bond)))
                     (setf (gethash side-chain-topology topology-to-residue) side-chain-residue)
                     (chem:add-plug side-chain-topology plug-name in-plug)
                     (push (cons (chem:get-name side-chain-residue) side-chain-topology) (gethash short-name side-chains))))
          (values core-topology side-chains topology-to-residue map-names-numbers))))))

(defun cartesian (l) 
  (cond ((null l) nil) 
        ((null (cdr l)) 
         (mapcar #'list (car l))) 
        (t (mapcan #'(lambda (x) (mapcar #'(lambda (y) (cons y x)) (car l))) 
                   (cartesian (cdr l))))))

(defvar *fep-time* 0)

(defun incf-fep-time ()
  (incf *fep-time*))

(defclass job-graph ()
  ((nodes :initform nil :accessor nodes)
   (edges :initform nil :accessor edges)))

(defclass fep-structure ()
  ((name :initarg :name :accessor name)
   (drawing :initarg :drawing :accessor drawing)
   (molecule :initarg :molecule :accessor molecule)
   (molecule-time :initform 0 :accessor molecule-time)
   (atom-order :initarg :atom-order :accessor atom-order)
   (atom-order-time :initform 0 :accessor atom-order-time)
   (am1-charges :initarg :am1-charges :accessor am1-charges)
   (am1-charges-time :initform 0 :accessor am1-charges-time)
   (am1-bcc-charges :initarg :am1-bcc-charges :accessor am1-bcc-charges)))

(defmethod print-object ((feps fep-structure) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (feps stream)
        (format stream "fep ~a" (string (name feps))))))

(defclass calculation ()
  ((receptors :initform nil :accessor receptors)
   (ligands :initarg :ligands :accessor ligands)
   (jobs :initarg :jobs :accessor jobs)))

(defun structure-name (tops)
  (let ((names (mapcar (lambda (x) (subseq (string (car x)) 1 3)) tops)))
    (intern (apply #'concatenate 'string (sort names #'string<)) :keyword)))

(defun set-coordinates-for-residue (res other-res)
  (chem:map-atoms
   nil
   (lambda (a)
     (let* ((other-atom (chem:atom-with-name other-res (chem:get-name a)))
            (pos (chem:get-position other-atom)))
       (chem:set-position a pos)))
   res))

(defun set-stereocenters (mol)
  (let ((stereocenters 
          (sort (cando:gather-stereocenters mol) #'string< :key #'chem:get-name)))
    (cando:set-stereoisomer-func stereocenters (constantly :S))))

(defun build-feps (core-topology side-topologys topology-to-residue map-names-numbers)
  (let ((top-groups nil))
    (maphash (lambda (name tops)
               (push tops top-groups))
             side-topologys)
    (let ((cartesian (cartesian top-groups)))
      (loop for cart in cartesian
            for name = (structure-name cart)
            for molecule = (chem:make-molecule name)
            for core-residue = (let ((res (chem:build-residue-single-name core-topology)))
                                 (set-coordinates-for-residue res (gethash core-topology topology-to-residue))
                                 res)
            collect (progn
                      (chem:add-matter molecule core-residue)
                      (loop for (name . top) in cart
                            for side-residue = (let ((res (chem:build-residue-single-name top)))
                                                 (set-coordinates-for-residue res (gethash top topology-to-residue))
                                                 res)
                            for in-plug = (chem:get-in-plug top)
                            for in-plug-name = (chem:get-name in-plug)
                            for out-plug-name = (chem:other-plug-name in-plug-name)
                            do (progn
                                 (chem:add-matter molecule side-residue)
                                 (maphash (lambda (name number)
                                            (let ((atom (chem:atom-with-name core-residue name)))
                                              #+(or)(format t "Atom: ~a   number: ~a~%" atom number)
                                              (chem:set-property atom :map number)))
                                          map-names-numbers)
                                 (chem:connect-residues core-topology
                                                        core-residue
                                                        out-plug-name
                                                        top
                                                        side-residue
                                                        in-plug-name)))
                      (set-stereocenters molecule)
                      (make-instance 'fep-structure :name name :drawing molecule))))))


(defun setup-ligands (feps sketch &key verbose)
  (let ((ligands (multiple-value-bind (core-topology side-chains topology-to-residue map-names-numbers)
                     (build-topologys sketch :verbose verbose)
                   (build-feps core-topology side-chains topology-to-residue map-names-numbers))))
    (setf (ligands feps) ligands)))

(defun layout-ligands (calculation &key (xdelta 15.0) (ydelta 15.0) (accessor 'drawing))
  (let* ((feps (copy-list (ligands calculation)))
         (numx (floor (sqrt (length feps))))
         (numy (1+ numx))
         (agg (chem:make-aggregate nil)))
    (loop for ix from 0 below numx
          for x = (* (- ix (floor (/ numx 2))) xdelta)
          do (loop for iy from 0 below numy
                   for y = (* (- iy (floor (/ numy 2))) ydelta)
                   for transform = (geom:make-m4-translate (geom:vec x y 0.0))
                   do (let ((fep (pop feps)))
                        (when fep
                          (let ((molecule (chem:matter-copy (funcall accessor fep))))
                            (chem:apply-transform-to-atoms molecule transform)
                            (chem:add-matter agg molecule))))))
    agg))

(defun layout-2 (name1 name2 feps &key (accessor 'drawing))
  (let ((mol1 (chem:matter-copy (funcall accessor (find name1 feps :key 'name))))
        (mol2 (chem:matter-copy (funcall accessor (find name2 feps :key 'name))))
        (agg (chem:make-aggregate 'two))
        (translate (geom:make-m4-translate (geom:vec 20.0 0.0 0.0))))
    (chem:add-matter agg mol2)
    (chem:apply-transform-to-atoms mol1 translate)
    (chem:add-matter agg mol1)
    agg))
    

(defun build-ligands (calculations)
  (let ((feps (ligands calculations)))
    (loop for fep in feps
          for mol = (chem:matter-copy (drawing fep))
          do (format t "Building ~a~%" (name fep))
          do (chem:fill-in-implicit-hydrogens mol)
          do (cando:build-unbuilt-hydrogens mol)
          do (funcall (find-symbol "ASSIGN-ATOM-TYPES" :leap) mol)
          do (setf (molecule fep) mol
                   (molecule-time fep) (incf-fep-time))
          do (format t "build-structures fep: ~a   map-atoms: ~a~%" fep (sorted-map-atoms fep)))))

(defun minimize-ligands (calculation)
  (let ((feps (ligands calculation)))
    (loop for fep in feps
          for mol = (molecule fep)
          for agg = (chem:make-aggregate nil)
          do (chem:add-matter agg mol)
          do (funcall (find-symbol "ASSIGN-ATOM-TYPES" :leap) agg)
          do (energy:minimize agg ))))



(defun sorted-map-atoms (fep)
  (let ((molecule (molecule fep))
        (map-atoms nil))
    (chem:map-atoms nil (lambda (a)
                          (when (chem:matter-get-property-or-default a :map nil)
                            (push (cons (chem:matter-get-property a :map) a) map-atoms)))
                    molecule)
    (sort map-atoms #'< :key #'car)))


(defun anchor-to-pose (moveable-atoms fixed-pose-atoms &key stereochemical-restraints)
  (mapc (lambda (fixed-pose-atom moveable-atom)
          (let* ((pos (chem:get-position fixed-pose-atom))
                 (anchor-restraint (core:make-cxx-object
                                    'chem:restraint-anchor
                                    :atom moveable-atom
                                    :position pos
                                    :weight 100.0)))
            (chem:clear-restraints moveable-atom)
            (chem:add-restraint moveable-atom anchor-restraint)
            ;; Assign stereochemical restraint if there is one that matches
            (let* ((map-id (chem:matter-get-property-or-default moveable-atom :map nil))
                   (stereochemistry-assignment (assoc map-id stereochemical-restraints)))
              (when stereochemistry-assignment
                (chem:set-stereochemistry-type moveable-atom :chiral)
                (format t "chem:get-stereochemistry-type ~a -> ~a~%" moveable-atom (chem:get-stereochemistry-type moveable-atom))
                (chem:set-configuration moveable-atom (cdr stereochemistry-assignment))
                (format t "Set stereochemistry of ~a to ~a~%" moveable-atom (cdr stereochemistry-assignment))))
            (format t "Anchored ~a to ~a restraints: ~a~%" moveable-atom pos (chem:all-restraints moveable-atom))))
        fixed-pose-atoms moveable-atoms))
                  
(defun pose-ligands (calculation fixed-atoms &key stereochemical-restraints)
  (let ((feps (ligands calculation)))
    (loop for fep in feps
          for molecule = (molecule fep)
          for moveable-atoms = (mapcar #'cdr (sorted-map-atoms fep))
          do (format t "fep ~a   moveable-atoms -> ~a~%" fep moveable-atoms)
          do (chem:superpose-one molecule moveable-atoms fixed-atoms)
          do (anchor-to-pose moveable-atoms fixed-atoms :stereochemical-restraints stereochemical-restraints)))
  (minimize-ligands calculation))



(defun build-initial-jobs (calculation &optional (connections 3))
  (let ((jobs (make-instance 'job-graph)))
    (let* ((unsorted-feps (copy-list (ligands calculation)))
           (sorted-feps (sort unsorted-feps #'< :key (lambda (x) (chem:number-of-atoms (drawing x))))))
      (loop for fep in sorted-feps
            for added = (nodes jobs)
            do (push fep (nodes jobs))
            do (loop for other in (subseq added 0 (min (length added) connections))
                     do (push (list fep other) (edges jobs)))))
    (setf (jobs calculation) jobs)))


(defun am1-file-name (fep type)
  (make-pathname :directory '(:relative "ligands") :name (format nil "am1-~a" (string (name fep))) :type type :defaults *default-pathname-defaults*))

(defun am1-calculation-complete (fep)
  "Return T if the AM1 charges for this FEP exist, are up to date and are available.
Otherwise return NIL."
  (let ((in-filename (am1-file-name fep "in"))
        (out-filename (am1-file-name fep "out"))
        (atom-order (atom-order fep)))
    (cond
      ((and (probe-file out-filename)
            (< (file-write-date in-filename) (file-write-date out-filename))
            (charges:read-am1-charges out-filename atom-order)))
      (t nil))))

(defun setup-am1-calculations (calculation &key (maxcyc 9999))
  (let ((feps (ligands calculation)))
    (let ((*default-pathname-defaults* (merge-pathnames "ligands/" *default-pathname-defaults*)))
      (ensure-directories-exist *default-pathname-defaults*)
      (format t "Creating am1 scripts in ~a~%" (truename *default-pathname-defaults*))
      (let ((all-target (make-string-output-stream))
            (targets (make-string-output-stream)))
        (loop for fep in feps
              for in-filename = (make-pathname :name (format nil "am1-~a" (string (name fep))) :type "in")
              for out-filename = (make-pathname :name (format nil "am1-~a" (string (name fep))) :type "out")
              for in-file = (open in-filename :direction :output)
              do (format targets "~a : ~a~%" (namestring out-filename) (namestring in-filename))
              do (format targets "~c$(AMBERHOME)/bin/sqm -O -i ~a -o ~a~%" #\tab (namestring in-filename) (namestring out-filename))
              do (format targets "~%")
              do (let ((atoms (charges:write-sqm-calculation in-file (molecule fep) :maxcyc maxcyc)))
                   (setf (atom-order fep) atoms))
              do (format all-target "~a " (namestring out-filename)))
        (with-open-file (makefile "makefile.charges" :direction :output)
          (format makefile "all : ~a~%" (get-output-stream-string all-target))
          (format makefile "~cecho Done all~%~%" #\tab)
          (format makefile "~a~%" (get-output-stream-string targets)))))))

(defun check-am1-calculations (calculation)
  (let* ((status-entries (loop for fep in (ligands calculation)
                               for name = (string (name fep))
                               for status = (am1-calculation-complete fep)
                               collect (list name status)))
         (sorted-entries (sort status-entries #'string< :key #'car)))
    (cl-jupyter-user:html (cl-markup:markup
                           (:table :border 3
                                   :cellpadding 4
                                   (loop for i below (length sorted-entries) by 5
                                         collect (cl-markup:markup
                                                  (:tr :align "right"
                                                       (loop for j from i below (+ i 5)
                                                             for entry = (elt sorted-entries j)
                                                             for name = (car entry)
                                                             for status = (cadr entry)
                                                             when entry
                                                               collect (cl-markup:markup
                                                                        (:td :bgcolor (if status "white" "red")
                                                                             (format nil "~a ~a" name (if status "Done" "Incomplete")))))))))))))
(defun read-am1-charges (calculation)
  (let ((*default-pathname-defaults* (merge-pathnames "ligands/" *default-pathname-defaults*))
        (count 0)
        (feps (ligands calculation)))
    (loop for fep in feps
       for out-filename = (make-pathname :name (format nil "am1-~a" (string (name fep))) :type "out")
       do (let ((am1-charge (charges::read-am1-charges out-filename (fep::atom-order (nth count feps)))))
            (setf (am1-charges fep) am1-charge))
       do (format t "file ~a done~%" (name fep))
       do (incf count))))


(defun calculate-am1-bcc-charges (calculation)
  (loop for fep in (ligands calculation)
        for am1-charges = (am1-charges fep)
        for bcc-corrections = (charges::calculate-bcc-corrections (fep::molecule fep))
        for am1-bcc-charges = (charges::calculate-am1-bcc-charges (am1-charges fep) bcc-corrections)
        do (setf (am1-bcc-charges fep) am1-bcc-charges)))


(defun complex (calculation ligand-name &optional (receptor-index 0))
  (let ((receptor (elt (receptors calculation) receptor-index))
        (ligand (find ligand-name (ligands calculation) :test #'eq :key #'name)))
    (let ((agg (chem:matter-copy receptor)))
      (cando:remove-molecules agg 'cando:ligand)
      (chem:add-molecule agg (fep:molecule ligand))
      agg)))


(defun relax-receptor-around-ligand (receptor-aggregate ligand-molecule &key (weight 100.0))
  "Relax the receptor around the ligand."
  (cando:do-atoms (atom receptor-aggregate)
    (let ((restraint-anchor (core:make-cxx-object
                             'chem:restraint-anchor
                             :atom atom
                             :position (chem:get-position atom)
                             :weight weight)))
      (chem:clear-restraints atom)
      (chem:add-restraint atom restraint-anchor)))
  (chem:add-matter receptor-aggregate ligand-molecule)
  (cando:do-atoms (atom ligand-molecule)
    (let ((restraint-anchor (core:make-cxx-object
                             'chem:restraint-anchor
                             :atom atom
                             :position (chem:get-position atom)
                             :weight weight)))
      (chem:clear-restraints atom)
      (chem:add-restraint atom restraint-anchor))))
