(in-package :asdf-user)

(defsystem "fep"
  :description "FEP setup code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on ( :cl-markup
                :cando
                :charges
                :leap
                :alexandria
                :cl-jupyter
                (:version :esrap "0.15")
                :parser.common-rules
                :PARSER.COMMON-RULES.OPERATORS
                :architecture.builder-protocol)
  :serial t
  :components
  ((:file "packages")
   (:file "jsme")
   (:file "graph2")
   (:file "fep")
   ))

