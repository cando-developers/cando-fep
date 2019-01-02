(in-package :ti)

(require :cl-ppcre)


(setf *default-pathname-defaults* #P"~/Dropbox/Shared Group Members/Shiho-Chris/working/fep-demo/fep/")

;;; What is this for????
(defparameter *vdw-bonded* "ifsc=1, scmask1=':1@H6', scmask2=':2@O1,H6'")


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

(defclass ti-file ()
  ((option :initarg :option :accessor option)
   (file :initarg :file :accessor file)))

(defun make-ti-file (option file)
  (make-instance 'ti-file :option option :file file))

(defclass step ()
  ((lam :initarg :lam :accessor lam)
   (step-info :initarg :step-info :accessor step-info)
   (amber-in :initarg :amber-in :accessor amber-in)
   (amber-in-file :initarg :amber-in-file :accessor amber-in-file)
   (inputs :initarg :inputs :accessor inputs)
   (outputs :initarg :outputs :accessor outputs)
   (makefile-clause :initarg :makefile-clause :accessor makefile-clause)))

(defclass step-info ()
  ((lam :initarg :lam :accessor lam)
   (unique-name :initform (gensym) :accessor unique-name)))

(defun make-step-info (lam)
  (make-instance 'step-info :lam lam))

(defclass heat-step (step) ())
(defclass ti-step (step) ())

(defun step-name (step)
  "This name is used to construct the directory for the step"
  (format nil "s-~,3f-~a" (lam step) (unique-name step)))

(defun make-heat-step (lam info)
  (make-instance 'heat-step
                 :lam lam
                 :step-info info
                 :inputs (list (make-ti-file "-i" "%esl%/heat.in")
                               (make-ti-file "-c" "%data%/ti-%e%.rst7")
                               (make-ti-file "-ref" "%data%/ti-%e%.rst7")
                               (make-ti-file "-p" "%data%/ti-%e%.parm7"))
                 :outputs (list (make-ti-file "-o" "%esl%/heat.out")
                                (make-ti-file "-inf" "%esl%/heat.info")
                                (make-ti-file "-e" "%esl%/heat.en")
                                (make-ti-file "-r" "%esl%/heat.rst7")
                                (make-ti-file "-x" "%esl%/heat.nc")
                                (make-ti-file "-l" "%esl%/heat.log"))
                 :amber-in *heat-in*
                 :amber-in-file "%esl%/heat.in"
                 :makefile-clause "%outputs% : %inputs%
	runcmd -- %inputs% -- %outputs% -- \\
	pmemd %option-inputs% \\
	  -O %option-outputs%"))

(defun make-ti-step (lam info)
  (make-instance 'ti-step
                 :lam lam
                 :step-info info
                 :inputs (list (make-ti-file "-i" "%esl%/ti.in")
                               (make-ti-file "-c" "%data%/heat.rst7")
                               (make-ti-file "-p" "%data%/ti-%e%.parm7"))
                 :outputs (list (make-ti-file "-o" "%esl%/ti001.out")
                                (make-ti-file "-inf" "%esl%/ti001.info")
                                (make-ti-file "-e" "%esl%/ti001.en")
                                (make-ti-file "-r" "%esl%/ti001.rst7")
                                (make-ti-file "-x" "%esl%/ti001.nc")
                                (make-ti-file "-l" "%esl%/ti001.log"))
                 :amber-in *ti-in*
                 :amber-in-file "%esl%/ti.in"
                 :makefile-clause "%outputs% : %inputs%
	runcmd -- %inputs% -- %outputs% -- \\
	pmemd %option-inputs% \\
	  -O %option-outputs%"))

(defclass ti-path ()
  ((lambdas :initform 11 :initarg :lambdas :accessor lambdas)
   (steps :initform nil :initarg :steps :accessor steps)
   (source-compound :initarg :source-compound :accessor source-compound)
   (target-compound :initarg :target-compound :accessor target-compound)))

(defun make-ti-path (lambdas start end)
  (make-instance 'ti-path :lambdas lambdas
                          :source-compound start
                          :target-compound end))

(defparameter *identical-lambda-delta* 0.001)
(defmethod maybe-add-step (ti-path step)
  (prog1
      (pushnew step (steps ti-path) :test (lambda (x y)
                                            (and (eq (class-of x) (class-of y))
                                                 (let ((delta (abs (- (lam x) (lam y)))))
                                                   (< delta *identical-lambda-delta*)))))
    (format t "Number of steps ~a~%" (length (steps ti-path)))))

(defun edge-string (ti-path)
  (format nil "~a-~a"
          (name (source-compound ti-path))
          (name (target-compound ti-path))))

(defun replace-all (dict in-script)
  (loop for script = in-script then script-result
        for (match . substitution) in dict
        do (format t "Replacing ~a with ~a~%" match substitution)
        do (setf script-result (cl-ppcre:regex-replace-all match script substitution))
        finally (return-from replace-all script-result)))

(defun generate-steps (ti-path)
  (let ((delta (/ 1.0 (1- (lambdas ti-path)))))
    (format t "delta  = ~a   lambdas = ~a~%" delta (lambdas ti-path))
    (loop for window from 0 below (lambdas ti-path)
          for lambda-val = 0.0 then (incf lambda-val delta)
          for step-info = (make-step-info lambda-val)
          for heat-step = (make-heat-step lambda-val step-info)
          for ti-step = (make-ti-step lambda-val step-info)
          do (format t "Adding step lambda ~f~%" lambda-val)
          do (maybe-add-step ti-path heat-step)
             (maybe-add-step ti-path ti-step))))
  
(defun generate-scripts (ti-path &key stage (directory *default-pathname-defaults*) makefile-stream)
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
            with edge-string = (edge-string ti-path)
            for esl = (format nil "~a/~a/~a" edge-string stage step-name)
            for substitutions = (list (cons "%e%" edge-string)
                                      (cons "%s%" stage)
                                      (cons "%l%" lambda-string)
                                      (cons "%lambdas,%" comma-separated-lambdas)
                                      (cons "%esl%" esl)
                                      (cons "%vdw-bonded%" *vdw-bonded*))
            for one-dir = (merge-pathnames (pathname esl) directory)
            for clause = (let* ((inputs (format nil "~{~a ~}" (mapcar (lambda (x) (file x)) (inputs step))))
                                (outputs (format nil "~{~a ~}" (mapcar (lambda (x) (file x)) (outputs step))))
                                (option-inputs (with-output-to-string (sout)
                                                 (mapc (lambda (x) (format sout "~a ~a " (option x) (file x))) (inputs step))))
                                (option-outputs (with-output-to-string (sout)
                                                  (mapc (lambda (x) (format sout "~a ~a " (option x) (file x))) (outputs step))))
                                (substitutions (list* (cons "%inputs%" inputs)
                                                      (cons "%outputs%" outputs)
                                                      (cons "%option-inputs%" option-inputs)
                                                      (cons "%option-outputs%" option-outputs)
                                                      substitutions)))
                           (replace-all substitutions (makefile-clause step)))
            for amber-in-file = (replace-all substitutions (amber-in-file step))
            for amber-in = (replace-all substitutions (amber-in step))
            do (format t "lambda-string = ~a  step = ~s~%" lambda-string step)
            do (ensure-directories-exist (merge-pathnames amber-in-file))
            do (with-open-file (fout (merge-pathnames amber-in-file) :direction :output :if-exists :supersede)
                 (format fout "~a~%" amber-in))
            do (format makefile-stream "~a~%~%" clause)))))


(defun generate-all-scripts (calculation)
  (let* ((*default-pathname-defaults* (merge-pathnames (pathname #p"meister/"))))
    (ensure-directories-exist *default-pathname-defaults*)
    (cando:save-cando calculation "calculation")
    (with-open-file (fout (make-pathname :name "makefile") :direction :output)
      (loop for edge in edges
            do (generate-scripts edge :makefile-stream fout :stage "ligand"))
      (loop for edge in edges
            do (generate-scripts edge :makefile-stream fout :stage "complex")))))

            
(defparameter *edges*
  (let* ((x1 (make-ti-compound "x1"))
         (x2 (make-ti-compound "x2"))
         (x3 (make-ti-compound "x3"))
         (path1 (make-ti-path 11 x1 x2))
         (_ (generate-steps path1))
         (path2 (make-ti-path 11 x2 x3))
         (_ (generate-steps path2))
         (calc (make-instance 'ti-calculation :compounds (list x1 x2 x3) :paths (list path1 path2))))
    calc))

*edges*

(defparameter *edges* (make-instance 'ti-calculation))
(cando:save-cando *edges* "/tmp/calc")

(cando:save-cando (make-instance 'foo) "/tmp/test")


(defclass bar () ()) (let ((*print-readably* t)) (print (make-instance 'bar)))

(defmethod print-object ((x foo) stream)
  (format stream "#S(foo)"))

(generate-all-scripts *edges*)

