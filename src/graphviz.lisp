(in-package :fepdot)

(defmethod draw-node (id (node fep:amber-job) stream)
  (format stream "~a [label = \"~a\" ];~%" id (fep:name (fep:script node))))

(defmethod draw-node (id (node fep:cpptraj-job) stream)
  (format stream "~a [label = \"~a\" ];~%" id (fep:name (fep:script node))))

(defmethod draw-node (id (node fep:argument) stream)
  (format stream "~a [label = \"~a/~a.~a\",shape=rectangle];~%" id (fep:morph-string (fep:morph (fep:node node))) (fep:name (fep:node node)) (fep:extension (fep:node node))))

(defmethod draw-edge (source-id target-id label stream)
  (format stream "~a -> ~a [label = \"~a\"];~%" source-id target-id label))

(defun maybe-draw-node (key node stream id-map)
  (unless key
    (error "key is nil for node ~s~%" node))
  (let ((id (gethash key id-map)))
    (if id
        id
        (let ((tid (gensym)))
          (draw-node tid node stream)
          (setf (gethash key id-map) tid)
          tid))))

(defun draw-one-job (job stream id-map)
  (let ((job-id (maybe-draw-node job job stream id-map)))
    (loop for input in (fep:inputs job)
          for drawn-input = (maybe-draw-node (fep:node input) input stream id-map)
          do (draw-edge drawn-input job-id (fep:option input) stream))
    (loop for output in (fep:outputs job)
          for drawn-output = (maybe-draw-node (fep:node output) output stream id-map)
          do (draw-edge job-id drawn-output (fep:option output) stream))))

(defun gather-jobs (job seen-jobs)
  (unless (gethash job seen-jobs)
    (setf (gethash job seen-jobs) t)
    (loop for output in (fep:outputs job)
          do (loop for user in (fep:users output)
                   do (gather-jobs user seen-jobs)))))

(defun draw-graph-stream (jobs stream)
  (let ((unique-jobs (make-hash-table))
        (id-map (make-hash-table)))
    (loop for job in jobs
          do (gather-jobs job unique-jobs))
    (format stream "digraph G {~%")
    (maphash (lambda (job dummy)
               (declare (ignore dummy))
               (draw-one-job job stream id-map))
             unique-jobs)
    (format stream "}~%")))
