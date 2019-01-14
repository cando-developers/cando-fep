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

(defun make-script-1-leap (calculation &key input-feps-file)
  (warn "in default-script-1-leap")
  (with-top-directory (calculation)
    (let (work-list)
      (loop for receptor in (receptors calculation)
            for jobs = (jobs calculation)
            ;; Write out the ligands
            do (format t "jobs -> ~s~%" jobs)
            do (loop for side in '(:ligand :complex)
                     for side-name = (format nil "~a-vdw-bonded" (string-downcase side))
                     do (loop for morph in (morphs jobs)
                              for source = (source morph)
                              for target = (target morph)
                              for parm-node = (make-instance 'morph-side-topology-file
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
                                                                                           :topology parm-node
                                                                                           :pdb pdb-node)
                                                                       :makefile-clause (standard-cando-makefile-clause script)))
                              do (multiple-value-bind (first-job last-job)
                                     (make-morph-side-prepare morph side
                                                              :input-topology-file parm-node
                                                              :input-coordinate-file coord-node)
                                   (let* ((press (output-file last-job :-r))
                                          (strip-job (make-morph-side-strip morph side
                                                                            :input-topology-file parm-node
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
                                                     :makefile-clause (standard-cando-makefile-clause script)))
                                     (push first-job work-list))))))
      work-list)))

(defmethod generate-jobs (calculation)
  (with-top-directory (calculation)
    (let* ((jupyter-job (make-instance 'jupyter-job))
           (am1-jobs (setup-am1-calculations jupyter-job calculation))
           (feps-precharge (make-instance 'feps-precharge-file)))
      (cando:save-cando calculation (node-pathname feps-precharge))
      (push (make-instance 'argument :option :feps-precharge :node feps-precharge) (outputs jupyter-job))
      (let* ((script (make-instance 'cando-script-file
                                    :name "charge"
                                    :script *cando-charge-script*))
             (feps-out (make-instance 'feps-postcharge-file)))
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
        (make-script-1-leap calculation :input-feps-file feps-out)
        ;; Do more preparation
        (generate-all-code calculation (list jupyter-job))
        jupyter-job))))
