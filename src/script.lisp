(in-package :fep)

(defclass fep-calculation (calculation)
  ((solvent-box :initform :tip3pbox :initarg :solvent-box :accessor solvent-box)
   (solvent-buffer :initform 12.0 :initarg :solvent-buffer :accessor solvent-buffer)
   (solvent-closeness :initform 0.75 :initarg :solvent-closeness :accessor solvent-closeness)
   (script-0-setup :initform 'default-script-0-setup :accessor script-0-setup)
   (script-1-leap :initform 'default-script-1-leap :accessor script-1-leap)))

(defmethod print-object ((obj fep-calculation) stream)
  (if *print-readably*
      (progn
        (format stream "#$(fep::fep-calculation ")
        (loop for slot in (clos:class-slots (find-class 'fep::fep-calculation))
              for slot-name = (clos:slot-definition-name slot)
              for initargs = (clos:slot-definition-initargs slot)
              if (and (car initargs) (slot-boundp obj slot-name))
                do (format stream "~s ~s " (car initargs) (slot-value obj slot-name)))
        (format stream ") "))
      (print-unreadable-object (obj stream)
        (format stream "fep-calculation"))))

#+(or)
(defun maybe-decharge-recharge-job (&key morph side press-job input-topology-file)
  (when (decharge-recharge morph)
    (let* ((press (output-file last-job :-r))
           (strip-job (make-morph-side-strip morph side
                                             :input-topology-file input-topology-file
                                             :input-coordinate-file press))
           (solvated (output-file strip-job :solvated))
           (source (output-file strip-job :source))
           (target (output-file strip-job :target))
           (script (make-instance 'morph-side-script
                                  :morph morph
                                  :side side
                                  :script *decharge-recharge-4-leap*
                                  :name "decharge-recharge"
                                  :extension "lisp"))
           (decharge-pdb (make-instance 'morph-side-pdb-file :morph morph :side side :name "decharge"))
           (decharge-topology (make-instance 'morph-side-topology-file :morph morph :side side :name "decharge"))
           (decharge-coordinates (make-instance 'morph-side-coordinate-file :morph morph :side side :name "decharge"))
           (recharge-pdb (make-instance 'morph-side-pdb-file :morph morph :side side :name "recharge"))
           (recharge-topology (make-instance 'morph-side-topology-file :morph morph :side side :name "recharge"))
           (recharge-coordinates (make-instance 'morph-side-coordinate-file :morph morph :side side :name "recharge")))
      (connect-graph
       (make-instance 'morph-side-cando-job
                      :morph morph
                      :side side
                      :inputs (arguments :feps input-feps-file
                                         :solvated solvated
                                         :source source
                                         :target target)
                      :outputs (arguments :decharge-pdb decharge-pdb
                                          :decharge-topology decharge-topology
                                          :decharge-coordinates decharge-coordinates
                                          :recharge-pdb recharge-pdb
                                          :recharge-topology recharge-topology
                                          :recharge-coordinates recharge-coordinates) 
                      :script script
                      :makefile-clause (standard-cando-makefile-clause script))))))


(defmacro powerloop ((&rest clauses) &rest final)
  (if (null clauses)
      `(progn ,@final)
    `(loop ,@(first clauses) do (powerloop (,@(rest clauses)) ,@final))))
 
(defun make-script-1-leap (calculation &key input-feps-file)
  (warn "in default-script-1-leap")
  (with-top-directory (calculation)
    (let (work-list morph-jobs)
      (powerloop
       ((for receptor in (receptors calculation)
             for jobs = (jobs calculation))
        (for morph in (morphs jobs)
             for source = (source morph)
             for target = (target morph)))
       (let (analysis-jobs)
         (loop
           for side in '(:ligand :complex)
           for side-name = (format nil "~a-vdw-bonded" (string-downcase side))
           for input-topology-file = (make-instance 'morph-side-topology-file
                                                    :morph morph :side side
                                                    :name side-name
                                                    :extension "parm7")
           for coord-node = (make-instance 'morph-side-coordinate-file
                                           :morph morph :side side
                                           :name side-name
                                           :extension "rst7")
           for pdb-node = (make-instance 'morph-side-pdb-file
                                         :morph morph :side side
                                         :name side-name)
           for script = (make-instance 'morph-side-script
                                       :morph morph
                                       :side side
                                       :script *solvate-addion-morph-side-script*
                                       :name "solvate-addion"
                                       :extension "lisp")
           for solvate-addion-job = (connect-graph
                                     (make-instance 'morph-side-cando-job
                                                    :morph morph
                                                    :side side
                                                    :script script
                                                    :inputs (arguments :input input-feps-file)
                                                    :outputs (arguments :coordinates coord-node
                                                                        :topology input-topology-file
                                                                        :pdb pdb-node)
                                                    :makefile-clause (standard-cando-makefile-clause script)))
           for morph-side-prepare-job = (make-morph-side-prepare morph side
                                                                 :input-topology-file input-topology-file
                                                                 :input-coordinate-file coord-node)
           for strip-job = (make-morph-side-strip morph side
                                                  :input-topology-file input-topology-file
                                                  :input-coordinate-file (output-file morph-side-prepare-job :-r))
           do (let (stage-jobs)
                (loop for stage in (case (stages morph)
                                     (1 '(:vdw-bonded))
                                     (3 '(:decharge :vdw-bonded :recharge))
                                     (otherwise (error "Illegal number of stages ~a - only 1 or 3 are allowed for morph ~s" (stages morph) morph)))
                      with source = (output-file strip-job :source)
                      with target = (output-file strip-job :target)
                      with solvated = (output-file strip-job :solvated)
                      for script-source = (ecase stage
                                            (:decharge *decharge*)
                                            (:vdw-bonded nil)
                                            (:recharge *recharge*))
                      for script = (unless (eq stage :vdw)
                                     (make-instance 'morph-side-stage-script-file :morph morph :side side :stage stage :script script-source
                                                                                  :name (string-downcase stage)))
                      for inputs = (ecase stage
                                     (:decharge (arguments :feps input-feps-file
                                                           :solvated solvated
                                                           :source source))
                                     (:vdw-bonded nil)
                                     (:recharge (arguments :feps input-feps-file
                                                           :solvated solvated
                                                           :target target)))
                      for outputs = (ecase stage
                                      (:decharge
                                       (arguments
                                        :decharge-pdb (make-instance 'morph-side-stage-pdb-file
                                                                     :morph morph
                                                                     :side side
                                                                     :stage stage
                                                                     :name "decharge")
                                        :decharge-topology (make-instance 'morph-side-stage-topology-file
                                                                          :morph morph
                                                                          :side side
                                                                          :stage stage
                                                                          :name "decharge")
                                        :decharge-coordinates (make-instance 'morph-side-stage-coordinates-file
                                                                             :morph morph
                                                                             :side side
                                                                             :stage stage
                                                                             :name "decharge")))
                                      (:vdw-bonded nil)
                                      (:recharge
                                       (arguments
                                        :recharge-pdb (make-instance 'morph-side-stage-pdb-file
                                                                     :morph morph
                                                                     :side side
                                                                     :stage stage
                                                                     :name "recharge")
                                        :recharge-topology (make-instance 'morph-side-stage-topology-file
                                                                          :morph morph
                                                                          :side side
                                                                          :stage stage
                                                                          :name "recharge")
                                        :recharge-coordinates (make-instance 'morph-side-stage-coordinates-file
                                                                             :morph morph
                                                                             :side side
                                                                             :stage stage
                                                                             :name "recharge"))))
                      for stage-job = (unless (eq stage :vdw-bonded)
                                        (connect-graph
                                         (make-instance 'morph-side-stage-cando-job
                                                        :morph morph :side side :stage stage :script script
                                                        :inputs inputs
                                                        :outputs outputs
                                                        :makefile-clause (standard-cando-makefile-clause script))))
                      do (let (decharge-jobs
                               vdw-jobs
                               recharge-jobs
                               (lambda-values (loop for window-index from 0 below (windows morph)
                                                    collect (/ (float window-index) (1- (windows morph))))))
                           (loop for lambda-value in lambda-values
                                 for lambda-label = (format nil "~5,3f" lambda-value)
                                 do (case stage
                                      (:decharge (let* ((input-topology-file (output-file stage-job :decharge-topology))
                                                        (input-coordinate-file (output-file stage-job :decharge-coordinates))
                                                        (heat-job (make-heat-ti-step morph side stage lambda-label lambda-values
                                                                                     :input-topology-file input-topology-file
                                                                                     :input-coordinate-file input-coordinate-file))
                                                        (input-coordinate-file (output-file heat-job :-r)))
                                                   (push (make-ti-step morph side stage lambda-label lambda-values
                                                                       :input-topology-file input-topology-file
                                                                       :input-coordinate-file input-coordinate-file)
                                                         decharge-jobs)))
                                      (:vdw-bonded (let* ((input-topology-file input-topology-file)
                                                          (input-coordinate-file (output-file morph-side-prepare-job :-r))
                                                          (heat-job (make-heat-ti-step morph side stage lambda-label lambda-values
                                                                                       :input-topology-file input-topology-file
                                                                                       :input-coordinate-file input-coordinate-file))
                                                          (input-coordinate-file (output-file heat-job :-r)))
                                                     (push (make-ti-step morph side stage lambda-label lambda-values
                                                                         :input-topology-file input-topology-file
                                                                         :input-coordinate-file input-coordinate-file)
                                                           vdw-jobs)))
                                      (:recharge (let* ((input-topology-file (output-file stage-job :recharge-topology))
                                                        (input-coordinate-file (output-file stage-job :recharge-coordinates))
                                                        (heat-job (make-heat-ti-step morph side stage lambda-label lambda-values
                                                                                     :input-topology-file input-topology-file
                                                                                     :input-coordinate-file input-coordinate-file))
                                                        (input-coordinate-file (output-file heat-job :-r)))
                                                   (push (make-ti-step morph side stage lambda-label lambda-values
                                                                       :input-topology-file input-topology-file
                                                                       :input-coordinate-file input-coordinate-file)
                                                         recharge-jobs)))))
                           (let ((script (make-instance 'morph-side-stage-script-file
                                                        :morph morph :side side :stage stage
                                                        :script *python-analyze*
                                                        :name "analyze"
                                                        :extension "py")))
                             (push (connect-graph
                                    (make-instance 'morph-side-stage-python-job
                                                   :morph morph :side side :stage stage :script script
                                                   :inputs (case stage
                                                             (:decharge (arguments :. (mapcar (lambda (job) (output-file job :-e)) decharge-jobs)))
                                                             (:vdw-bonded (arguments :. (mapcar (lambda (job) (output-file job :-e)) vdw-jobs)))
                                                             (:recharge (arguments :. (mapcar (lambda (job) (output-file job :-e)) recharge-jobs))))
                                                   :outputs (arguments :stage-analysis (make-instance 'morph-side-stage-file
                                                                                                      :morph morph :side side :stage stage
                                                                                                      :name "dvdl" :extension "dat"))
                                                   :makefile-clause (standard-makefile-clause "analyze.py -- %inputs%")))
                                   stage-jobs))))
                ;; combine stage-jobs
                (let ((side-script (make-instance 'morph-side-script
                                                  :morph morph :side side
                                                  :script *combine-stages*
                                                  :name "side-script"
                                                  :extension "lisp")))
                  (push (connect-graph
                         (make-instance 'morph-side-cando-job
                                        :morph morph :side side
                                        :inputs (arguments :. (mapcar (lambda (job) (output-file job :stage-analysis)) stage-jobs))
                                        :outputs (arguments :side-analysis (make-instance 'morph-side-file
                                                                                          :morph morph :side side
                                                                                          :name "side-analysis" :extension "dat"))
                                        :makefile-clause (standard-cando-makefile-clause side-script)))
                        analysis-jobs))))
         ;; combine analysis-jobs
         (let ((script (make-instance 'morph-script
                                      :morph morph
                                      :script *combine-sides*
                                      :name "morph-script"
                                      :extension "lisp")))
           (push (connect-graph
                  (make-instance 'morph-cando-job
                                 :morph morph
                                 :inputs (arguments :. (mapcar (lambda (job) (output-file job :side-analysis)) analysis-jobs))
                                 :outputs (arguments :morph-analysis (make-instance 'morph-file
                                                                                    :morph morph
                                                                                    :name "morph-analysis" :extension "dat"))
                                 :makefile-clause (standard-cando-makefile-clause script)))
                 morph-jobs))))
      morph-jobs)))

(defmethod generate-jobs (calculation)
  (let ((*default-pathname-defaults* (merge-pathnames (top-directory calculation) *default-pathname-defaults*)))
    (format t "FUCK ~s~%" *default-pathname-defaults*)
    (let* ((jupyter-job (make-instance 'jupyter-job))
           (am1-jobs (setup-am1-calculations jupyter-job calculation))
           (feps-precharge (make-instance 'feps-file :name "precharge")))
      (fep:save-feps calculation (node-pathname feps-precharge))
      (push (make-instance 'argument :option :feps-precharge :node feps-precharge) (outputs jupyter-job))
      (let* ((script (make-instance 'cando-script-file
                                    :name "charge"
                                    :script *cando-charge-script*))
             (feps-out (make-instance 'feps-file :name "postcharge")))
        (connect-graph
         (make-instance 'cando-job
                        :inputs (apply #'arguments
                                       :feps feps-precharge
                                       (loop for am1-job in am1-jobs
                                             for output = (output-file am1-job :-o)
                                             append (list :input output)))
                        :outputs (arguments :output feps-out)
                        :script script
                        :makefile-clause (standard-cando-makefile-clause script)))
        (let ((morph-jobs (make-script-1-leap calculation :input-feps-file feps-out)))
          ;; Do more preparation
          (generate-all-code calculation (list jupyter-job) (mapcar (lambda (job) (output-file job :morph-analysis)) morph-jobs)))
        jupyter-job))))
