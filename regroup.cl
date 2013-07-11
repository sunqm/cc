;;; -*- Mode: Common Lisp; -*-
;;;
;;; 

(load "utility.cl")

; regroup polynomial
(defun sum? (poly)
  (eql '+ (car poly)))
(defun prod? (poly)
  (eql '* (car poly)))
(defun remove-prod-operator (prod)
  (if (eql '* (car prod))
      (cdr prod)
      prod))
(defun intersectionp (&rest lsts)
  (let ((lst1 (car lsts))
        (lst-rest (cdr lsts)))
    (some (lambda (item)
            (every (lambda (lst) (member item lst :test #'equal))
                   lst-rest))
          lst1)))

;;; When the second most common term does not overlap with the most
;;; common term, it will be found in the next search.  This fn tries
;;; to find out the common terms which are mutually exclusive.
;;; if no common-term is found, return nil
(defun find-best-common-term (n-bests lst-prod fn-term-id)
  (labels ((count-obj (obj)
             (count-if (lambda (prod) (member (funcall fn-term-id obj)
                                              prod :test #'equal))
                       (mapcaar fn-term-id lst-prod)))
           ;; find unique terms which are reasonable for regrouping.
           (find-uniq-term (filter)
             (remove-duplicates (apply #'append (mapcar filter lst-prod))
                                :key fn-term-id :test #'equal))
           (keep-prod-which-has (keys)
             (lambda (prod)
               (let ((keys-id (mapcar fn-term-id keys))
                     (prod-id (mapcar fn-term-id prod)))
                 (if (subsetp keys-id prod-id :test #'equal)
                     (remove-if (lambda (x) (member (funcall fn-term-id x) keys
                                                    :key fn-term-id :test #'equal))
                                prod)))))
           (build-match-list (obj count matches)
             (if (or (eql count 1) (null obj))
                 matches ; no common term is found, do not add to the match lists
                 (cons (list count obj) matches)))
           (find-most (n)
             (if (eql n 1)
                 (multiple-value-bind (obj count)
                     (most #'count-obj (find-uniq-term #'identity))
                   (build-match-list obj count '()))
                 (let ((matches (find-most (1- n))))
                   (if (null matches)
                       '()
                       (let* ((keys (mapcar #'second matches))
                              (objs (find-uniq-term (keep-prod-which-has keys))))
                         (if objs
                             (multiple-value-bind (obj count) (most #'count-obj objs)
                               (build-match-list obj count matches))
                             matches)))))))
    (if (null lst-prod)
        '()
        (reverse (find-most n-bests)))))

(defun find-most-common-labels (n lst-prod)
  (find-best-common-term n (mapcar #'remove-prod-operator lst-prod) #'identity))

;;; return two sets: the first set includes obj and the second set collects the rest
(defun subgroup-prods (term lst-prod)
  (flet ((contain-obj (prod)
           (member term prod :test #'equal))
         (single-term? (prod)
           (last-one? (remove-prod-operator prod))))
    (let* ((set1 (remove-if-not #'contain-obj lst-prod))
           (set-rest (remove-if #'contain-obj lst-prod))
           (sub-set1 (mapcar (lambda (prod) (remove term prod :test #'equal))
                             set1))
           (singles (mapcar #'second (remove-if-not #'single-term? sub-set1)))
           (multi-s (remove-if #'single-term? sub-set1)))
      (values singles multi-s set-rest))))

;;; given a list of products to be sumed up, regroup them to
;;; a list of polynomials; each polynomial is sth like (+ (* a b (+ c ..)))
;;; 2 is enough for most case. Very low possibility that 3 most common
;;; terms are adjoined in the same term-products.
(defparameter *common-term-search-trial* 2)
;;; find all possible regroup candidates
(defun regroup-iter (poly)
  (cond ((null poly) (list '()))
        ((last-one? poly) (list poly))
        (t (flet ((expand (term)
                    (multiple-value-bind (singles multi-s set-rest)
                        (subgroup-prods term poly)
                      (flet ((comb (regp-mults regp-rest)
                               `((* ,term (+ ,@singles ,@regp-mults))
                                 ,@regp-rest)))
                        (mapnest #'comb
                                 (regroup-iter multi-s)
                                 (regroup-iter set-rest))))))
             (let ((commons (find-most-common-labels
                             *common-term-search-trial* poly)))
               (if (null commons)
                   (list poly)
                   (mapcan (lambda (common-term) (expand (cadr common-term)))
                           commons)))))))
(defun regroup (poly)
  (mapcar (lambda (res) (cons '+ res))
          (regroup-iter poly)))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; define the filter function which considers memory usage
;;; It restricts the number of open lines
(defparameter *max-ext-lines* 4)
; fixme: control the usage of memory
(defun filter-for-tight-mem (ampprod)
  (flet ((ext-int-line (amp)
           (let* ((n-int (count-if #'internal-line?
                                   (apply #'append (content-of amp))))
                  (n-ext (- (* 2 (num-nodes-in-amp amp)) n-int)))
             (values n-ext n-int))))
    (let ((tot-ext (apply #'+ (mapcar #'ext-int-line ampprod))))
      (remove-if (lambda (amp)
                   (multiple-value-bind (n-ext n-int) (ext-int-line amp)
                     ;(print (list tot-ext (+ (- tot-ext n-ext) n-int) amp))
                     (> (+ (- tot-ext n-ext) n-int) *max-ext-lines*)))
                 ampprod))))

(defun amp-id-with-e1<=>e2 (amp)
  (let ((sid1 (id-of-amp amp))
        (sid2 (id-of-amp (replace-amp-index
                          #'line-indexed? '((1 . 2) (2 . 1))
                          amp))))
    (if (list<= sid1 sid2)
        sid1
        sid2)))
(defun find-most-common-amp (n ampprod-lst)
    (find-best-common-term n (mapcar #'filter-for-tight-mem ampprod-lst)
                         #'amp-id-with-e1<=>e2))

(defun toggle-amp-line-int<->ext (amp line)
  (cons (tag-of amp)
        (mapcaar (lambda (l)
                   (if (equal l line)
                       (toggle-line-int<->ext l)
                       l))
                   (content-of amp))))
(defun cut-op-from-amp (op amp)
  (reduce (lambda (amp line)
            (if (internal-line? line)
                (let ((conn-line (toggle-line-excite<->dexcite line)))
                  (toggle-amp-line-int<->ext amp conn-line))
                amp))
          (apply #'append (content-of op))
          :initial-value amp))
;;; if amp is a member of ampprod, remove amp and toggle the
;;; internal/external line in ampprod
(defun cut-amp-from-ampprod (amp ampprod)
  (let ((subprod (remove (id-of-amp amp) ampprod
                         :key #'id-of-amp :test #'equal :count 1)))
    (mapcar (lambda (x) (cut-op-from-amp amp x))
            subprod)))

(defun subgroup-ampprods (term ampprod-lst)
  (flet ((label-term-via (fn-id)
           (lambda (prod)
             (member (funcall fn-id term) prod
                     :key fn-id :test #'equal)))
         (single-term? (prod)
           (last-one? prod)))
    (let* ((set1 (remove-if-not (label-term-via #'amp-id-with-e1<=>e2) ampprod-lst))
           (set-rest (remove-if (label-term-via #'amp-id-with-e1<=>e2) ampprod-lst))
           (sub1 (remove-if-not (label-term-via #'id-of-amp) set1))
           (sub2 (mapcar #'swap-ampprod-e1<->e2
                         (remove-if (label-term-via #'id-of-amp) set1)))
           (cut-set1 (mapcar (lambda (ampprod)
                               (cut-amp-from-ampprod term ampprod))
                             (append sub2 sub1)))
           (singles (mapcar #'car (remove-if-not #'single-term? cut-set1)))
           (multi-s (remove-if #'single-term? cut-set1)))
      (values singles multi-s set-rest))))

(defun regroup-amp-iter (poly)
  (cond ((null poly) (list '()))
        ((last-one? poly)
         (list `(,(cons '* (car poly)))))
        (t (flet ((expand (term)
                    (multiple-value-bind (singles multi-s set-rest)
                        (subgroup-ampprods term poly)
                      (flet ((comb (regp-mults regp-rest)
                               `((* ,term (+ ,@singles ,@regp-mults))
                                 ,@regp-rest)))
                        (mapnest #'comb
                                 (regroup-amp-iter multi-s)
                                 (regroup-amp-iter set-rest))))))
             (let ((commons (find-most-common-amp
                             *common-term-search-trial* poly)))
               (if (null commons)
                   (list (mapcar (lambda (x) (cons '* x)) poly))
                   (mapcan (lambda (common-term) (expand (cadr common-term)))
                           commons)))))))








(defun count-fp (no nv))
