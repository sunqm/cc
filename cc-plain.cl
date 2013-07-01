;;; -*- Mode: Common Lisp; -*-
;;;
;;; Description: based on antisymmetrized Goldstone diagram

(load "utility.cl")

(defun attach-tag (tag x)
  (cons tag x))
(proclaim '(inline attach-tag))
(defmacro tag-of (x)
  `(car ,x))
(defmacro content-of (x)
  `(cdr ,x))

;;; type:
;;; he+   ->-|
;;; pe+   -<-|
;;; he-   |->-
;;; pe-   |-<-
;;; hi+   contracted ->-|
;;; pi+   contracted -<-|
;;; hi-   contracted |->-
;;; pi-   contracted |-<-
;;;
;;; (t_{ijk}^{abc} i a j b k c)

(defun op-copy (n op)
  (loop repeat n
       collect op))

; amp == (t op1 op2 ...)
; n-excitation amplitude without linking
(defun make-new-amp (n)
  (attach-tag 't (flatten (op-copy n '(he+ pe+)))))
; t1 = (t he+ pe+)
; t2 = (t he+ pe+ he+ pe+)

;;; contraction of one hole and an amplitude
(defun contract-hole-amp (amp)
  (let ((ops (replace-oncep (lambda (op) (if (eql op 'he+) 'hi+))
                            (content-of amp))))
    (if ops
        (attach-tag 't ops))))
(defun contract-particle-amp (amp)
  (let ((ops (replace-oncep (lambda (op) (if (eql op 'pe+) 'pi+))
                            (content-of amp))))
    (if ops
        (attach-tag 't ops))))
;;; try to avoid redundant contraction
(defun contract-hole-amp* (amp)
  (let* ((first-ctr t)
         (ops (replace-oncep (lambda (op)
                               (cond ((eql op 'pi+) (setf first-ctr nil))
                                     ((and first-ctr (eql op 'pe+)) 'pi+)))
                             (content-of amp))))
    (if ops
        (attach-tag 't ops))))
(defun contract-particle-amp* (amp)
  (let ((ops (replace-oncep (lambda (op) (if (eql op 'pe+) 'pi+))
                            (content-of amp))))
    (if ops
        (attach-tag 't ops))))

(defun contract-hole-ampprod (ampprod)
  (let ((last-amp '()))
    (mapreplace (lambda (amp)
                  (unless (equal last-amp amp)
                    (setf last-amp amp)
                    (contract-hole-amp amp)))
                ampprod)))
(defun contract-particle-ampprod (ampprod)
  (let ((last-amp '()))
    (mapreplace (lambda (amp)
                  (unless (equal last-amp amp)
                    (setf last-amp amp)
                    (contract-particle-amp amp)))
                ampprod)))

(defun contract-ops-ampprods (ops ampprod-lst)
  (flet ((hole-func (ampprod-lst op)
           (if (eql op 'he-)
               (let ((res (mapcan #'contract-hole-ampprod ampprod-lst)))
                 (if res res ampprod-lst))
               ampprod-lst))
         (particle-func (ampprod-lst op)
           (if (eql op 'pe-)
               (let ((res (mapcan #'contract-particle-ampprod ampprod-lst)))
                 (if res res ampprod-lst))
               ampprod-lst)))
    (reduce #'particle-func ops ; first reduce all holes to avoid sorting
            :initial-value (reduce #'hole-func ops :initial-value ampprod-lst))))

(defun id-amp (amp &optional (oplist '(hi- pi- hi+ pi+ he- pe- he+ pe+)))
  (let ((shift (length oplist)))
    (reduce (lambda (v x) (+ x (* shift v)))
            (mapcar (lambda (x) (position x oplist))
                    (sort (content-of amp)
                          #'< :key (lambda (x) (position x oplist))))
            :initial-value 0)))
(defun sort-ampprod (ampprod)
  (sort ampprod
        (lambda (x y)
          (let ((len-x (length x))
                (len-y (length y)))
            (if (eq len-x len-y)
                (< (id-amp x) (id-amp y))
                (< len-x len-y))))))

;;; result list of contract-ops-ampprods may need to be sorted
(defun contract-ops-ampprods-uniq (ops ampprod-lst)
  (remove-duplicates (contract-ops-ampprods ops ampprod-lst)
                     :test #'equal))

(defun connected-amp? (amp)
  (some (lambda (op) (member op '(hi- pi- hi+ pi+)))
        (content-of amp)))
(defun connected-ampprod? (ampprod)
  (every #'connected-amp? ampprod))

(defun ext-line-to-int-line (op)
  (case op
    (he- 'hi-)
    (pe- 'pi-)
    (he+ 'hi+)
    (pe+ 'pi+)
    (otherwise op)))
(defun contract-h2e-ampprod (h2e ampprod)
  (let ((ctr-h2e (maptree (lambda (op)
                            (case op (he- 'hi-)
                                     (pe- 'pi-)
                                     (otherwise op)))
                          h2e))
        (ctr-amps (contract-ops-ampprods-uniq (content-of h2e)
                                              (list ampprod))))
    (mapcar (lambda (amps) (cons ctr-h2e amps))
            (remove-if-not #'connected-ampprod? ctr-amps))))

(defparameter *h2e-ops*
  '((g he- he+ he+ pe+)
    (g pe- pe+ he+ pe+)
    (g he- he+ he- he+)
    (g pe- pe+ pe- pe+)
    (g he- he+ pe- pe+)
    (g he- pe- he- he+)
    (g he- pe- pe- pe+)
    (g he- pe- he- pe-)))    
(defparameter *h1e-h2e*
  (append '((f he- he+)
            (f pe- pe+)
            (f he- pe-)
            (f he+ pe+))
          *h2e-ops*))

(defun gen-amps-list (n-lst)
  (if (null n-lst)
      '()
      (let* ((tnew (make-new-amp (car n-lst)))
             (t1 (list tnew))
             (t2 (cons tnew t1))
             (t3 (cons tnew t2))
             (t4 (cons tnew t3))
             (amps-lst (gen-amps-list (cdr n-lst)))
             (new-lst (mapcan (lambda (tx)
                                (mapcar (lambda (amps) (append tx amps))
                                        amps-lst))
                              (list t1 t2 t3)))
             (mod-lst (remove-if (lambda (ts) (> (length ts) 4))
                                 new-lst)))
        (append (list t1 t2 t3 t4) mod-lst))))

(defun count-tot-lines (pred ampprod)
  (apply #'+ (mapcar (lambda (amp) (count-if pred (content-of amp)))
                     ampprod)))
(defun count-excite-lines (ampprod)
  (count-tot-lines (lambda (x) (member x '(he+ pe+)))
                   ampprod))
(defun count-dexcite-lines (ampprod)
  (count-tot-lines (lambda (x) (member x '(he- pe-)))
                   ampprod))

(defun gen-diagrams-w/o-index (n-lst)
  (flet ((d-ext-lines (ops ampprod)
           (- (count-excite-lines (cons ops ampprod))
              (count-dexcite-lines (list ops)))))
    (let ((avail-exts (mapcar (lambda (x) (+ x x)) n-lst)))
      (mapcan (lambda (ampprod)
                (mapcan (lambda (ops)
                          (if (member (d-ext-lines ops ampprod) avail-exts)
                              (contract-h2e-ampprod ops ampprod)))
                        *h1e-h2e*))
              (gen-amps-list n-lst)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; adding indices
;;; e.g. assign 2 to internal hole line (hi+ . 2)
;;;

;;; rline is an indexed operator == (op index)
(defun make-rline (op idx)
  (cons op idx))
(defun symb-of-rline (op)
  (car op))
(defun idx-of-rline (op)
  (cdr op))

; return nil if not replaced
(defun label-amp (test idx amp)
  (let ((ops (replace-oncep (lambda (op) (if (funcall test op)
                                             (make-rline op idx)))
                            (content-of amp))))
    (if ops
        (attach-tag (tag-of amp) ops))))
(defun label-amp-contract-hole (idx amp)
  (label-amp (lambda (op) (or (eql op 'hi-) (eql op 'hi+))) idx amp))
(defun label-amp-contract-particle (idx amp)
  (label-amp (lambda (op) (or (eql op 'pi-) (eql op 'pi+))) idx amp))

;;; if not contracted, return the input ampprod
(defun label-contract-pair (label-func idx h-ampprod)
  (cons (funcall label-func idx (car h-ampprod))
        (replace-once (lambda (amp) (funcall label-func idx amp))
                      (cdr h-ampprod))))

;;; return a product, sum over the indices which appear twice
(defun label-lines (ampprod)
  (let ((reg 0))
    (flet ((label-ctr (ts op)
             (case op
               (hi- (label-contract-pair #'label-amp-contract-hole
                                         (incf reg) ts))
               (pi- (label-contract-pair #'label-amp-contract-particle
                                         (incf reg) ts))
               (otherwise ts)))
           (label-ex (symb)
             (lambda (amp)
               (let ((pairs (mapcar (lambda (op)
                                      (if (eql op symb)
                                          (make-rline op (incf reg))
                                          op))
                                    (content-of amp))))
                 (attach-tag (tag-of amp) pairs)))))
      (mapcar (label-ex 'pe+)
              (mapcar (label-ex 'he+)
                      (reduce #'label-ctr (content-of (car ampprod))
                              :initial-value ampprod))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; interperate the diagrams
;;;

;;; equivalent internal lines give a factor of 1/2
;;; return the factor 1, 1/2 or 1/4
(defun factor-int-line-pair (h-ampprod)
  (flet ((count-eq-line (symb)
           (apply #'* (mapcar (lambda (amp)
                                (if (eql 2 (count symb (content-of amp) :key #'car))
                                    .5 1))
                              (cdr h-ampprod)))))
    (* (count-eq-line 'hi+)
       (count-eq-line 'pi+))))

;;; equivalent external lines must connect to the same vertex
(defun find-rline-vertex (rline h-ampprod)
  (find-if (lambda (amp) (member rline amp :test #'equal))
           h-ampprod))
(defun ext-line-eql? (ex1 ex2 h-ampprod)
  ;;; same symbol and same vertex
  (and (eql (symb-of-rline ex1) (symb-of-rline ex2))
       (equal (find-rline-vertex ex1 h-ampprod)
              (find-rline-vertex ex2 h-ampprod))))

(defun vertex-equiv? (amp1 amp2)
  (flet ((get-symb (x)
           (if (atom x) x (symb-of-rline x))))
    (equal (mapcar #'get-symb (content-of amp1))
           (mapcar #'get-symb (content-of amp2)))))
(defun ext-line-symm-eql? (ex1 ex2 h-ampprod)
  ;;; same symbol on equivalent vertex
  (vertex-equiv? (find-rline-vertex ex1 h-ampprod)
                 (find-rline-vertex ex2 h-ampprod)))

(defun factor-of-symm-vertices (h-ampprod)
  (let ((ampprod (cdr h-ampprod)))
    (case (length ampprod)
      (2 (if (vertex-equiv? (first ampprod) (second ampprod)) .5 1))
      (3 ;todo
       (and (equal (first ampprod) (second ampprod))
              (equal symm1?)))
      (4 ;todo
       (and (equal (first ampprod) (second ampprod))
            (equal (third ampprod) (fourth ampprod))))
      (otherwise 1))))

;;; return the connected rlines
(defun track-rline (rline h-ampprod)
  (let ((finds (list rline)))
    (labels ((find-friend (rline) ; find the rline on the same node
               (let* ((v (content-of (find-rline-vertex rline h-ampprod)))
                      (pos (position rline v :test #'equal)))
                 (if (evenp pos)
                     (nth (1+ pos) v)
                     (nth (1- pos) v))))
             (conn-int-line (rline)
               (let ((s (symb-of-rline rline))
                     (i (idx-of-rline rline)))
                 (case s
                   (hi- (make-rline 'hi+ i))
                   (hi+ (make-rline 'hi- i))
                   (pi- (make-rline 'pi+ i))
                   (pi+ (make-rline 'pi- i))
                   (otherwise nil))))
             (searching (rline)
               (let ((next-rline (find-friend rline)))
                 (setf finds (cons next-rline finds))
                 (if (member (symb-of-rline next-rline)
                             '(hi- hi+ pi- pi+))
                     (let ((nnext (conn-int-line next-rline)))
                       (setf finds (cons nnext finds))
                       (searching nnext))
                     finds))))
      (searching rline))))

(defun count-loops (ampprod)
  ;;; fixme: + inner loops
  (count-tot-lines (lambda (rline)
                     (eql (symb-of-rline rline) 'he+))
                   ampprod))
(defun count-hole-lines (ampprod)
  (count-tot-lines (lambda (rline)
                     (member (symb-of-rline rline) '(he- he+ hi- hi+)))
                   ampprod))
(defun hole-loop-sign (ampprod)
  (if (evenp (+ (count-hole-lines ampprod)
                (count-loops ampprod)))
      1
      -1))

(defun collect-ext-holes (ampprod)
  (remove-if-not #'identity
                 (mapcar (lambda (amp)
                           (find-if (lambda (rline)
                                      (eql (symb-of-rline rline) 'he+))
                                    (content-of amp)))
                         ampprod)))

;;; elem-lst has at least one item
(defun permutation-sets (elem-lst &key (test #'equal))
  (if (last-one? elem-lst)
      (list elem-lst)
      (mapcan (lambda (e)
                (mapcar (lambda (s) (cons e s))
                        (permutation-sets
                         (remove e elem-lst :count 1 :test #'equal))))
              (remove-duplicates elem-lst :test test :from-end t))))
;;; even permutation and odd permutation
(defun permutation-even-odd-sets (elem-lst &key (test #'equal))
  (labels ((per-it (elems)
             (if (last-one? elems)
                 (list (list elems) '())
                 (reduce (lambda (eo-sets e)
                           (let* ((eset (first eo-sets))
                                  (oset (second eo-sets))
                                  (p (per-it (remove e elems :count 1 :test #'equal)))
                                  (s1 (mapcar (lambda (s) (cons e s)) (first p)))
                                  (s2 (mapcar (lambda (s) (cons e s)) (second p))))
                             (if (evenp (position e elems))
                                 (list (append s1 eset) (append s2 oset))
                                 (list (append s2 eset) (append s1 oset)))))
                         (remove-duplicates elems :test test :from-end t)
                         :initial-value '(() ())))))
    (mapcar #'reverse (per-it elem-lst))))

;;; symmetric vertices will cancel against the permutation of inequivalent ext-lines
(defun permutation-ext-line (h-ampprod)
  (let* ((h-ext (collect-ext-holes h-ampprod))
         (p-ext (mapcar #'car (track-rline h-ext h-ampprod))))
    (flet ((ext-eql (x1 x2)
             (or (ext-line-eql? x1 x2 h-ampprod)
                 (ext-line-symm-eql? x1 x2 h-ampprod))))
      (list (permutation-even-odd-sets h-ext :test #'ext-eql)
            (permutation-even-odd-sets p-ext :test #'ext-eql)))))
