(in-package :fep)

(defclass fep-calculation (calculation)
  ((solvent-box :initform :tip3pbox :initarg :solvent-box :accessor solvent-box)
   (solvent-buffer :initform 12.0 :initarg :solvent-buffer :accessor solvent-buffer)
   (solvent-closeness :initform 0.75 :initarg :solvent-closeness :accessor solvent-closeness)
   (script-0-setup :initform 'default-script-0-setup :accessor script-0-setup)
   (script-1-leap :initform 'default-script-1-leap :accessor script-1-leap)))

(defun default-script-0-setup (calculation)
  (leap:source "leaprc.water.tip3p")
  (leap:load-off "solvents.lib")
  (leap:load-off "atomic_ions.lib"))
  
(defun default-script-1-leap (calculation)
  (warn "in default-script-1-leap")
  (let (work-list)
    (loop for receptor in (receptors calculation)
          for jobs = (jobs calculation)
          ;; Write out the ligands
          do (format t "jobs -> ~s~%" jobs)
          do (loop for morph in (morphs jobs)
                   for source = (source morph)
                   for target = (target morph)
                   for source-molecule = (chem:matter-copy (molecule source))
                   for target-molecule = (chem:matter-copy (molecule target))
                   for ligands = (cando:combine source-molecule target-molecule)
                   for ligands-name = (format nil "ligands-vdw-bonded-~a-~a"
                                              (string (name source))
                                              (string (name target)))
                   for parm-name = (make-instance 'morph-side-topology-file
                                                  :morph morph :side ':ligand
                                                  :name "ligands-vdw-bonded"
                                                  :extension "parm7")
                   for coord-name = (make-instance 'morph-side-coordinate-file
                                                   :morph morph :side ':ligand
                                                   :name "ligands-vdw-bonded"
                                                   :extension "rst7")
                   do (setf (ligands-name morph) ligands-name)
                      (leap:solvate-box ligands
                                        (leap.core:lookup-variable (solvent-box calculation))
                                        (solvent-buffer calculation)
                                        (solvent-closeness calculation))
                      (leap.add-ions:add-ions ligands :|Cl-| 0)
                      (chem:save-pdb ligands (target-pathname :pdb calculation ligands-name))
                      (ensure-directories-exist (node-pathname parm-name))
                      (ensure-directories-exist (node-pathname coord-name))
                      (leap.topology:save-amber-parm-format ligands
                                                            (node-pathname parm-name)
                                                            (node-pathname coord-name))
                      (multiple-value-bind (first-job last-job)
                          (make-morph-side-prepare morph :ligand
                                                   :input-topology-file parm-name
                                                   :input-coordinate-file coord-name)
                        (let ((press (output-file last-job :-r)))
                          (make-morph-side-strip morph :ligand
                                                 :input-topology-file parm-name
                                                 :input-coordinate-file press))
                        (push first-job work-list)))
             ;; Write out the complexes
          do (loop for morph in (morphs jobs)
                   for source = (source morph)
                   for target = (target morph)
                   for source-molecule = (molecule source)
                   for target-molecule = (molecule target)
                   for ligands = (cando:combine (chem:matter-copy source-molecule)
                                                (chem:matter-copy target-molecule))
                   for complex-name = (format nil "complex-~a-vdw-bonded-~a-~a"
                                              (string (chem:get-name receptor))
                                              (string (name source))
                                              (string (name target)))
                   for complex = (cando:combine (chem:matter-copy receptor)
                                                ligands)
                   for parm-name = (make-instance 'morph-side-topology-file
                                                  :morph morph :side ':complex
                                                  :name "complex-vdw-bonded"
                                                  :extension "parm7")
                   for coord-name = (make-instance 'morph-side-coordinate-file
                                                   :morph morph :side ':complex
                                                   :name "complex-vdw-bonded"
                                                   :extension "rst7")
                   do (setf (complex-name morph) complex-name)
                      (leap:solvate-box complex
                                        (leap.core:lookup-variable (solvent-box calculation))
                                        (solvent-buffer calculation)
                                        (solvent-closeness calculation))
                      (leap.add-ions:add-ions complex :|Cl-| 0)
                      (chem:save-pdb complex (target-pathname :pdb calculation complex-name))
                      (ensure-directories-exist (node-pathname parm-name))
                      (ensure-directories-exist (node-pathname coord-name))
                      (leap.topology:save-amber-parm-format complex
                                                            (node-pathname parm-name)
                                                            (node-pathname coord-name))
                      (multiple-value-bind (first-job last-job)
                          (make-morph-side-prepare morph :complex
                                                   :input-topology-file parm-name
                                                   :input-coordinate-file coord-name)
                        (push first-job work-list)
                        (let ((press (output-file last-job :-r)))
                          (make-morph-side-strip morph :complex
                                                 :input-topology-file parm-name
                                                 :input-coordinate-file press))
                        )))
    work-list))

(defmethod generate-files (calculation)
  (funcall (script-0-setup calculation) calculation)
  (funcall (script-1-leap calculation) calculation)
  ;; Do more preparation
  ;;;(generate-all-scripts calculation)
  )
