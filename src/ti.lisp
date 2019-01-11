(in-package :fep)

;;; What is this for????
(defparameter *vdw-bonded* "ifsc=1, scmask1=':1@H6', scmask2=':2@O1,H6'")

(defparameter *prepare-min-in*
"  minimisation
 &cntrl
   imin = 1, ntmin = 2,
   maxcyc = 100,
   ntpr = 20, ntwe = 20,
   ntb = 1,
   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, ifsc = 1, clambda = 0.5, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = '%timask1%', timask2 = '%timask2%',
   scmask1 = '%scmask1%', scmask2 = '%scmask2%'
 /
 &ewald
 / 
")

(defparameter *prepare-heat-in*
  "heating
 &cntrl
   imin = 0, nstlim = 10000, irest = 0, ntx = 1, dt = 0.002,
   nmropt = 1,
   ntt = 1, temp0 = 300.0, tempi = 5.0, tautp = 1.0,
   ntb = 1,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, ifsc = 1, clambda = 0.5, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = '%timask1%', timask2 = '%timask2%',
   scmask1 = '%scmask1%', scmask2 = '%scmask2%'
 /
 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 8000,                                      
   value1 = 5.0, value2 = 300.0
 /

 &wt type = 'END'
 /

")

(defparameter *prepare-press-in*
  "pressurising
 &cntrl
   imin = 0, nstlim = 10000, irest = 1, ntx = 5, dt = 0.002,
   ntt = 1, temp0 = 300.0, tautp = 1.0,
   ntp = 1, pres0 = 1.0, taup = 2.0,
   ntb = 2,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, ifsc = 1, clambda = 0.5, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = '%timask1%', timask2 = '%timask2%',
   scmask1 = '%scmask1%', scmask2 = '%scmask2%'
 /
 &ewald
 / 

")


(defparameter *cpptraj-strip*
"trajin %side%_prepare/press.rst7

# remove the two ligands and keep the rest
strip \":1,2\"
outtraj %side%_solvated.pdb onlyframes 1

# extract the first ligand
unstrip
strip \":2-999999\"
outtraj %side%_bnz.pdb onlyframes 1

# extract the second ligand
unstrip
strip \":1,3-999999\"
outtraj %side%_phn.pdb onlyframes 1
")


(defparameter *heat-in* 
  "heating
 &cntrl
   imin = 0, nstlim = 10000, irest = 0, ntx = 1, dt = 0.002,
   ntt = 1, temp0 = 300.0, tempi = 50.0, tautp = 1.0,
   ntc = 2, ntf = 1,
   ntb = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   nmropt = 1,
   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, clambda = %l%, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = ':1', timask2 = ':2',
   %vdw-bonded%
 /

 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 8000,                                      
   value1 = 50.0, value2 = 300.0
 /

 &wt type = 'END'
 /

")

(defparameter *ti-in*
  "TI simulation
 &cntrl
   imin = 0, nstlim = 100000, irest = 1, ntx = 5, dt = 0.002,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   ntc = 2, ntf = 1,
   ntb = 2,
   ntp = 1, pres0 = 1.0, taup = 2.0,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 10000, ntpr = 10000, ntwr = 20000,

   icfe = 1, clambda = %l%, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 1,
   ifmbar = 1, bar_intervall = 1000, mbar_states = 11,
   mbar_lambda = %lambdas,%
   timask1 = ':1', timask2 = ':2',
   %vdw-bonded%
 /

 &ewald
 / 

")

(defclass ti-calculation ()
  ((compounds :initarg :compounds :accessor compounds)
   (paths :initarg :paths :accessor paths)))

(defclass ti-compound ()
  ((name :initarg :name :accessor name)))

(defun make-ti-compound (name)
  (make-instance 'ti-compound :name name))

;;; ------------------------------------------------------------
;;;
;;; File node types
;;;
;;; side means either :ligand or :complex
;;; lambda means the lambda value from 0.0 to 1.0
;;; morph refers to the morph between two compounds of a compound graph

(defclass node-file ()
  ((definers :initarg :definers :initform nil :accessor definers)
   (users :initarg :users :initform nil :accessor users)
   (directory% :initarg :directory% :accessor directory%)
   (name :initarg :name :accessor name)
   (extension :initarg :extension :accessor extension)))


(defgeneric substitutions (calculation node-file))

(defclass morph-side-file (node-file)
  ((morph :initarg :morph :accessor morph)
   (side :initarg :side :accessor side))
  (:documentation "Files that only depend on morph and side"))


(defclass morph-side-unknown-file (morph-side-file)
  ()
  (:documentation "Files with no specific purpose that depend on the morph/side"))

(defclass morph-side-pdb-file (morph-side-file)
  ()
  (:documentation "PDB files"))

(defclass morph-side-script (morph-side-file)
  ((script :initarg :script :accessor script))
  (:default-initargs
   :extension "in")
  (:documentation "This is the input script file"))

(defmethod substitutions (calculation (node-file morph-side-script))
  (let ((morph (morph node-file)))
    (list (cons "%e%" (morph-string morph))
          (cons "%s%" (string-downcase (side node-file)))
          (cons "%timask1%" (calculate-timask calculation (source morph)))
          (cons "%timask2%" (calculate-timask calculation (target morph)))
          (cons "%scmask1%" (calculate-scmask calculation (source morph)))
          (cons "%scmask2%" (calculate-scmask calculation (target morph))))))

(defclass morph-side-restart-file (morph-side-file)
  ()
  (:default-initargs
   :extension "rst7")
  (:documentation "Restart files depending only on morph/side"))

(defclass morph-side-coordinate-file (morph-side-file)
  ()
  (:default-initargs
   :extension "crd")
  (:documentation "Coordinate files depending only on morph/side"))

(defclass morph-side-topology-file (morph-side-file)
  ()
  (:documentation "Topology files depending only on morph/side"))

(defclass morph-side-trajectory-file (morph-side-file)
  ()
  (:documentation "Trajectory files depending only on morph/side"))

(defclass morph-side-lambda-file (morph-side-file)
  ((lambda% :initarg :lambda% :accessor lambda%))
  (:documentation "A type of file node identified by a side (ligand vs complex),
a lambda value and an morph." ))

(defun morph-side-lambda-unknown-file (morph-side-lambda-file)
  ()
  (:documentation "A file that has an unknown purpose - use this as a placeholder until we know what
its for and then create a new class for it."))

(defclass morph-side-lambda-script (ti-morph-side-lambda-file)
  ((script :initarg :script :accessor script))
  (:default-initargs
   :extension "in")
  (:documentation "This is the input script file"))

(defclass morph-side-lambda-trajectory-file (morph-side-lambda-file)
  ()
  (:documentation "This is a trajectory file"))

(defclass morph-side-lambda-trajectory-file (morph-side-lambda-file)
  ()
  (:documentation "This is a restart file"))


(defun make-morph-side-lambda-script (&key name script)
  (make-instance 'ti-morph-side-lambda-script :name name :script script))

(defclass ti-morph-coordinates-file (ti-morph-file)
  ()
  (:documentation "This is a coordinates file"))

(defclass ti-morph-topology-file (ti-morph-file)
  ()
  (:documentation "This is a topology file"))

(defgeneric node-pathname (node))

(defmethod node-pathname ((node node-file))
  (make-pathname :name (name node)
                 :type (extension node)))

(defmethod node-pathname ((node morph-side-file))
  (merge-pathnames (call-next-method)
                   (make-pathname :directory (list
                                              :relative
                                              (string (morph-string (morph node)))
                                              (string-downcase (string (side node)))))))

(defun option-extension (option)
  (case option
    (:-i "in")
    (:-o "out")
    (:-inf "info")
    (:-e "en")
    (:-r "rst7")
    (:-x "nc")
    (:-l "log")
    (otherwise (error "Unknown sander/pmemd option ~a" option))))


;;; ------------------------------------------------------------
;;;
;;;




;;; ------------------------------------------------------------
;;;
;;;

(defclass job ()
  ((script :initarg :script :accessor script)
   (inputs :initform nil :initarg :inputs :accessor inputs)
   (outputs :initform nil :initarg :outputs :accessor outputs)
   (makefile-clause :initarg :makefile-clause :accessor makefile-clause)))

(defclass amber-job ()
  ())

(defclass cpptraj-job ()
  ())

(defclass morph-side-job (job)
  ((morph :initarg :morph :accessor morph)
   (side :initarg :side :accessor side))
  (:documentation "An AMBER job that only depends on the morph and side"))

(defclass morph-side-amber-job (morph-side-job amber-job)
  ())

(defclass morph-side-cpptraj-job (morph-side-job cpptraj-job)
  ())

(defclass ti-step (amber-job)
  ((step-info :initarg :step-info :accessor step-info)))

(defun output-file (amber-job option)
  (loop for output-arg in (outputs amber-job)
        when (eq option (option output-arg))
          do (return-from output-file (node output-arg)))
  (error "Could not find option ~a in ~a" option (outputs amber-job)))

(defclass step-info ()
  ((unique-step-name :initform (gensym) :accessor unique-step-name)))

(defun make-step-info (lam)
  (make-instance 'step-info :lam lam))

(defclass heat-step (step) ())

(defclass ti-step (step)
  ((lam :initarg :lam :accessor lam)))
(defclass ti-step-info (step-info)
  ((lam :initarg :lam :accessor lam)))

(defun make-ti-step-info (lam)
  (make-instance 'ti-step-info :lam lam))


(defun step-name (step)
  "This name is used to construct the directory for the step"
  (format nil "s-~,3f-~a" (lam step) (unique-step-name step)))

(defun make-morph-side-lambda-file (name extension)
  (make-ti-file :pathname (make-pathname :name (string name) :type (string extension))))

(defun connect-graph (amber-job)
  (let ((script (script amber-job))
        (inputs (inputs amber-job))
        (outputs (outputs amber-job)))
    (push amber-job (users script))
    (loop for input in inputs
          do (format t "input -> ~s~%" input)
             (format t "(users input) -> ~s~%" (users input))
          do (push amber-job (users (node input))))
    (loop for output in outputs
          do (format t "output -> ~s~%" output)
             (format t "(users output) -> ~s~%" (users output))
          do (push amber-job (definers (node output))))
    amber-job))

(defclass argument ()
  ((option :initarg :option :accessor option)
   (node :initarg :node :accessor node)))

(defmethod print-object ((obj argument) stream)
  (print-unreadable-object (obj stream)
    (format stream "~a ~a" (option obj) (node obj))))

(defmethod users ((arg argument))
  (users (node arg)))

(defmethod definers ((arg argument))
  (definers (node arg)))

(defun arguments (&rest args)
  "Later we will add sanity checks to the args - for now just accumulate a list"
  (loop for cur = args then (cddr cur)
        for option = (first cur)
        for node = (second cur)
        until (null cur)
        do (format t "option: ~a node: ~a~%" option node)
        while option
        collect (make-instance 'argument :option option :node node)))


(defun make-morph-side-prepare-step (morph side &key name script input-coordinate-file input-topology-file)
  (unless input-coordinate-file (error "You must provide an input-coordinate-file"))
  (unless input-topology-file (error "You must provide an input-topology-file"))
  (connect-graph
   (make-instance 'morph-side-amber-job
                  :morph morph
                  :side side
                  :script (make-instance 'morph-side-script :morph morph :side side :name name :script script)
                  :inputs (arguments :-c input-coordinate-file ; min.rst7
                                     :-ref input-coordinate-file ; ../ligands_vdw_bonded.rst7
                                     :-p input-topology-file) ; $prmtop
                  :outputs (arguments :-o (make-instance 'morph-side-unknown-file
                                                         :morph morph :side side
                                                         :name name :extension "out") ;  "%epl%/heat.out")
                                      :-inf (make-instance 'morph-side-unknown-file
                                                           :morph morph :side side
                                                           :name name :extension "info") ; (make-ti-file "-inf" "%epl%/heat.info")
                                      :-e (make-instance 'morph-side-unknown-file
                                                         :morph morph :side side
                                                         :name name :extension "en") ; (make-ti-file "-e" "%epl%/heat.en")
                                      :-r (make-instance 'morph-side-restart-file
                                                         :morph morph :side side
                                                         :name name) ; (make-ti-file "-r" "%epl%/heat.rst7")
                                      :-x (make-instance 'morph-side-trajectory-file
                                                         :morph morph :side side
                                                         :name name :extension "nc") ;  (make-ti-file "-x" "%epl%/heat.nc")
                                      :-l (make-instance 'morph-side-unknown-file
                                                         :morph morph :side side
                                                         :name name :extension "log") ; (make-ti-file "-l" "%epl%/heat.log"))
                                      )
                  :makefile-clause "%outputs% : %inputs%
	runcmd -- %inputs% -- %outputs% -- \\
	pmemd %option-inputs% \\
	  -O %option-outputs%")))

(defun make-morph-side-prepare (morph side &rest args &key input-coordinate-file input-topology-file)
  (unless input-coordinate-file (error "You must provide an input-coordinate-file"))
  (unless input-topology-file (error "You must provide an input-topology-file"))
  (let* ((min (make-morph-side-prepare-step morph side
                                            :name "min"
                                            :script *prepare-min-in*
                                            :input-coordinate-file input-coordinate-file
                                            :input-topology-file input-topology-file))
         (min.rst (output-file min :-r))
         (heat (make-morph-side-prepare-step morph side
                                             :name "heat"
                                             :script *prepare-heat-in*
                                             :input-coordinate-file min.rst
                                             :input-topology-file input-topology-file))
         (heat.rst (output-file heat :-r))
         (press (make-morph-side-prepare-step morph side
                                              :name "press"
                                              :script *prepare-press-in*
                                              :input-coordinate-file heat.rst
                                              :input-topology-file input-topology-file)))
    (values min press)))


(defun make-morph-side-strip (morph side &key input-coordinate-file input-topology-file)
  (format t "make-morph-side-strip ~a ~a ~a ~a~%" morph side input-coordinate-file input-topology-file)
  (make-instance 'morph-side-cpptraj-job
                 :morph morph
                 :side side
                 :script (make-instance 'morph-side-script :morph morph :side side :name "strip" :script *cpptraj-strip*)
                 :inputs (arguments :-p input-topology-file :coordinates input-coordinate-file)
                 :outputs (arguments :solvated (make-instance 'morph-side-pdb-file
                                                              :morph morph
                                                              :side side
                                                              :name "solvated"
                                                              :extension "pdb")
                                     :source (make-instance 'morph-side-pdb-file
                                                            :morph morph
                                                            :side side
                                                            :name "source"
                                                            :extension "pdb")
                                     :target (make-instance 'morph-side-pdb-file
                                                            :morph morph
                                                            :side side
                                                            :name "target"
                                                            :extension "pdb"))))

(defun make-heat-ti-step (lam info &key input-coordinate-file input-topology-file)
  (make-instance 'heat-step
                 :lam lam
                 :step-info info
                 :script (make-morph-side-lambda-script :name "heat" :script *heat-in*) ; "%epl%/heat.in"
                 :inputs (arguments :-c input-coordinate-file ; (make-ti-file "-c" "%data%/ti-%e%.rst7")
                                    :-ref input-coordinate-file ; (make-ti-file "-ref" "%data%/ti-%e%.rst7")
                                    :-p input-topology-file) ; (make-ti-file "-p" "%data%/ti-%e%.parm7")
                 :outputs (arguments :-o (make-instance 'morph-side-lambda-unknown-file :name "heat" :extension "out") ;  "%epl%/heat.out")
                                     :-inf (make-instance 'morph-side-lambda-unknown-file :name "heat" :extension "info") ; (make-ti-file "-inf" "%epl%/heat.info")
                                     :-e (make-instance 'morph-side-lambda-unknown-file :name "heat" :extension "en") ; (make-ti-file "-e" "%epl%/heat.en")
                                     :-r (make-instance 'morph-side-lambda-restart-file :name "heat") ; (make-ti-file "-r" "%epl%/heat.rst7")
                                     :-x (make-instance 'morph-side-lambda-trajectory-file :name "heat" :extension "nc") ;  (make-ti-file "-x" "%epl%/heat.nc")
                                     :-l (make-instance 'morph-side-lambda-unknown-file :name "heat" :extension "log") ; (make-ti-file "-l" "%epl%/heat.log"))
                                     )
                 :makefile-clause "%outputs% : %inputs%
	runcmd -- %inputs% -- %outputs% -- \\
	pmemd %option-inputs% \\
	  -O %option-outputs%"))

(defun make-ti-step (lam info &key input-coordinate-file input-topology-file)
  (let ((result (make-instance 'ti-step
                               :lam lam
                               :step-info info
                               :script (make-instance 'morph-side-lambda-script :name "ti" :script *ti-in*)
                               :inputs (arguments :-c input-coordinate-file
                                                  :-p input-topology-file)
                               :outputs (arguments :-o (make-instance 'morph-side-lambda-unknown-file :name "ti001" :extension "out") ; list (make-ti-file "-o" "%epl%/ti001.out")
                                                   :-inf (make-instance 'morph-side-lambda-unknown-file :name "ti001" :extension "info") ; (make-ti-file "-inf" "%epl%/ti001.info")
                                                   :-e (make-instance 'morph-side-lambda-unknown-file :name "ti001" :extension "en") ; (make-ti-file "-e" "%epl%/ti001.en")
                                                   :-r (make-instance 'morph-side-lambda-restart-file :name "ti001") ; (make-ti-file "-r" "%epl%/ti001.rst7")
                                                   :-x (make-instance 'morph-side-lambda-trajectory-file :name "ti001" :extension "nc") ; (make-ti-file "-x" "%epl%/ti001.nc")
                                                   :-l (make-instance 'morph-side-lambda-unknown-file :name "ti001" :extension "log") ; (make-ti-file "-l" "%epl%/ti001.log"))
                                                   :makefile-clause "%outputs% : %inputs%
	runcmd -- %inputs% -- %outputs% -- \\
	pmemd %option-inputs% \\
	  -O %option-outputs%"))))
    result))

(defclass ti-path ()
  ((lambdas :initform 11 :initarg :lambdas :accessor lambdas)
   (steps :initform nil :initarg :steps :accessor steps)
   (source-compound :initarg :source-compound :accessor source-compound)
   (target-compound :initarg :target-compound :accessor target-compound)))

(defparameter *identical-lambda-delta* 0.001)
(defmethod maybe-add-step (ti-path step)
  (prog1
      (pushnew step (steps ti-path) :test (lambda (x y)
                                            (and (eq (class-of x) (class-of y))
                                                 (let ((delta (abs (- (lam x) (lam y)))))
                                                   (< delta *identical-lambda-delta*)))))
    (format t "Number of steps ~a~%" (length (steps ti-path)))))

(defun replace-all (dict in-script)
  (loop for script = in-script then script-result
        for (match . substitution) in dict
        do (setf script-result (cl-ppcre:regex-replace-all match script substitution))
        finally (return-from replace-all script-result)))

(defun generate-steps (ti-path)
  (let ((delta (/ 1.0 (1- (lambdas ti-path)))))
    (format t "delta  = ~a   lambdas = ~a~%" delta (lambdas ti-path))
    (loop for window from 0 below (lambdas ti-path)
          for lambda-val = 0.0 then (incf lambda-val delta)
          for step-info = (make-ti-step-info lambda-val)
          for heat-step = (make-heat-step lambda-val step-info)
          for ti-step = (make-ti-step lambda-val step-info)
          do (format t "Adding step lambda ~f~%" lambda-val)
          do (maybe-add-step ti-path heat-step)
             (maybe-add-step ti-path ti-step))))

(defun make-ti-path (lambdas start end)
  (let ((path (make-instance 'ti-path :lambdas lambdas
                                      :source-compound start
                                      :target-compound end)))
    (generate-steps path)
    path))


(defun makefile-substitutions (calculation job)
  (let (inputs option-inputs outputs option-outputs)
    (loop for input in (inputs job)
          do (pushnew (namestring (node-pathname (node input))) inputs :test #'string=)
             (when (char= (char (string (option input)) 0) #\-)
               (push (string-downcase (option input)) option-inputs)
               (push (namestring (node-pathname (node input))) option-inputs)))
    (loop for output in (outputs job)
          do (pushnew (namestring (node-pathname (node output))) outputs :test #'string=)
             (when (char= (char (string (option output)) 0) #\-)
               (push (string-downcase (option output)) option-outputs)
               (push (namestring (node-pathname (node output))) option-outputs)))
    (list (cons "%inputs%" (format nil "~{~a ~}" (reverse inputs)))
          (cons "%outputs%" (format nil "~{~a ~}" (reverse outputs)))
          (cons "%option-inputs%" (format nil "~{~a ~}" (reverse option-inputs)))
          (cons "%option-outputs%" (format nil "~{~a ~}" (reverse option-outputs))))))
    

(defmethod generate-code (calculation (job amber-job) makefile visited-nodes)
  ;; Generate script
  (let* ((script (script job))
         (raw-script (script script))
         (substituted-script (replace-all (substitutions calculation script) raw-script)))
    (with-open-file (fout (ensure-directories-exist (node-pathname script)) :direction :output :if-exists :supersede)
      (write-string substituted-script fout))
    (let* ((raw-makefile-clause (makefile-clause job))
           (makefile-substitutions (makefile-substitutions calculation job))
           (substituted-makefile-clause (replace-all makefile-substitutions raw-makefile-clause)))
      (write-string substituted-makefile-clause makefile)
      (terpri makefile)
      (terpri makefile)))
  (loop for output in (outputs job)
        do (loop for child in (users output)
                 unless (gethash child visited-nodes)
                   do (setf (gethash child visited-nodes) t)
                      (generate-code calculation child makefile visited-nodes))))



(defun generate-all-code (calculation work-list)
  (let ((*default-pathname-defaults* (merge-pathnames #P"jobs/")))
    (let ((visited-nodes (make-hash-table))
          (makefile-pathname (ensure-directories-exist (merge-pathnames "makefile"))))
      (format t "Writing makefile to ~a~%" (translate-logical-pathname makefile-pathname))
      (with-open-file (makefile makefile-pathname :direction :output :if-exists :supersede)
        (loop for job in work-list
              do (generate-code calculation job  makefile visited-nodes))))))

#+(or)
(defun generate-scripts (ti-path &key side (directory *default-pathname-defaults*) makefile-stream)
  (let* ((steps (steps ti-path))
         (sorted-steps (sort (copy-list steps) #'< :key (lambda (x) (lam x))))
         (unique-lambdas (let (infos)
                           (loop for s in sorted-steps
                                 for info = (step-info s)
                                 do (pushnew info infos :test #'eq))
                           (mapcar #'lam (nreverse infos))))
         (comma-separated-lambdas (format nil "~{~,3f,~} " unique-lambdas)))
    (with-output-to-string (s)
      (loop for step in sorted-steps
            for lambda-val = (lam step)
            for lambda-string = (format nil "~,3f" lambda-val)
            for step-info = (step-info step)
            for step-name = (step-name step-info)
            with morph-string = (morph-string ti-path)
            for epl = (format nil "~a/~a/~a" morph-string side step-name)
            for substitutions = (list (cons "%e%" morph-string)
                                      (cons "%s%" side)
                                      (cons "%l%" lambda-string)
                                      (cons "%lambdas,%" comma-separated-lambdas)
                                      (cons "%epl%" epl)
                                      (cons "%vdw-bonded%" *vdw-bonded*))
            for one-dir = (merge-pathnames (pathname epl) directory)
            for clause = (let* ((inputs (format nil "~{~a ~}" (mapcar (lambda (x) (file x)) (inputs step))))
                                (outputs (format nil "~{~a ~}" (mapcar (lambda (x) (file x)) (outputs step))))
                                (option-inputs (with-output-to-string (sout)
                                                 (mapc (lambda (x) (format sout "~a ~a "
                                                                           (option x)
                                                                           (file x))) (inputs step))))
                                (option-outputs (with-output-to-string (sout)
                                                  (mapc (lambda (x) (format sout "~a ~a " (option x) (file x))) (outputs step))))
                                (substitutions (list* (cons "%inputs%" inputs)
                                                      (cons "%outputs%" outputs)
                                                      (cons "%option-inputs%" option-inputs)
                                                      (cons "%option-outputs%" option-outputs)
                                                      substitutions)))
                           (replace-all substitutions (makefile-clause step)))
            for script = (replace-all substitutions (script step))
            for amber-in = (replace-all substitutions (amber-in step))
            do (format t "lambda-string = ~a  step = ~s~%" lambda-string step)
               (ensure-directories-exist (merge-pathnames script))
               (with-open-file (fout (merge-pathnames script) :direction :output :if-exists :supersede)
                 (format fout "~a~%" amber-in))
               (format makefile-stream "~a~%~%" clause)))))



  
#+(or)
(defun generate-all-scripts (calculation work-list)
  (let* ((*default-pathname-defaults* (merge-pathnames (pathname #p"meister/"))))
    (ensure-directories-exist *default-pathname-defaults*)
    ;; (cando:save-cando calculation "calculation")
    (loop for job in work-list
          for script = (script job)
          for script = (specialize-script (script script))
          do (with-open-file (fout (node-pathname script) :direction :output :if-exists :supersede)
               (format fout "~a~%" script)))
    #+(or)
    (with-open-file (fout (make-pathname :name "makefile") :direction :output)
      (loop for morph in (morphs (jobs calculation))
            for morph-path = (make-ti-path (ti-lambdas calculation) (source morph) (target morph))
            do (generate-scripts morph :makefile-stream fout :side "ligand"))
      (loop for morph in (morphs (jobs calculation))
            do (generate-scripts morph :makefile-stream fout :side "complex")))))



#|
(defparameter *morphs*
  (let* ((x1 (make-ti-compound "x1"))
         (x2 (make-ti-compound "x2"))
         (x3 (make-ti-compound "x3"))
         (path1 (make-ti-path 11 x1 x2))
         (_ (generate-steps path1))
         (path2 (make-ti-path 11 x2 x3))
         (_ (generate-steps path2))
         (calc (make-instance 'ti-calculation :compounds (list x1 x2 x3) :paths (list path1 path2))))
    calc))

*morphs*

(defparameter *morphs* (make-instance 'ti-calculation))
(cando:save-cando *morphs* "/tmp/calc")

(cando:save-cando (make-instance 'foo) "/tmp/test")


(defclass bar () ()) (let ((*print-readably* t)) (print (make-instance 'bar)))

(defmethod print-object ((x foo) stream)
  (format stream "#S(foo)"))

(generate-all-scripts *morphs*)

|#
