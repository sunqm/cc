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
           (term-id-eql (term1 term2)
             (equal (funcall fn-term-id term1)
                    (funcall fn-term-id term2)))
           ;; find unique terms which are reasonable for regrouping.
           (find-uniq-term (filter)
             (remove-duplicates (apply #'append (mapcar filter lst-prod))
                                :test #'term-id-eql))
           (keep-prod-with (keys)
             (lambda (prod)
               (let ((keys-id (mapcar fn-term-id keys))
                     (prod-id (mapcar fn-term-id prod)))
                 (if (subsetp keys-id prod-id)
                     (remove-if (lambda (x) (member x keys :test #'term-id-eql))
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
                              (objs (find-uniq-term (keep-prod-with keys))))
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
;;; define the filter function which considers memory usage
;;; It restricts the number of open lines
(defparameter *avial-ext-lines* 4)
; fixme: control the usage of memory
(defun filter-for-tight-mem (ampprod)
  (flet ((ext-int-line (amp)
           (let* ((n-int (count-if (lambda (x)
                                     (if (listp x)
                                         (member (symb-of x) '(hi+ pi+ hi- pi-))))
                                   (apply #'append (content-of amp))))
                  (n-ext (- (* 2 (length (content-of amp))) n-int)))
             (values n-ext n-int))))
    (let ((tot-ext (apply #'+ (mapcar #'ext-int-line ampprod))))
      (remove-if (lambda (amp)
                   (multiple-value-bind (n-ext n-int) (ext-int-line amp)
                     (> (+ (- tot-ext n-ext) n-int) *avial-ext-lines*)))
                 ampprod))))
(defun amp-type-iden? (amp1 amp2)
  (let ((tid1 (symb-id-of-amp amp1))
        (tid2 (symb-id-of-amp amp2)))
    (equal tid1 tid2)))

(defun find-most-common-amp (n ampprod-lst)
  (find-best-common-term n (mapcar #'filter-for-tight-mem ampprod-lst)
                         #'symb-id-of-amp))

(defun subgroup-ampprods (term ampprod-lst)
  (flet ((contain-obj (prod)
           (member term prod :test #'amp-type-iden?))
         (single-term? (prod)
           (last-one? prod)))
    ;;; fixme: change line-symb
    (let* ((set1 (remove-if-not #'contain-obj ampprod-lst))
           (set-rest (remove-if #'contain-obj ampprod-lst))
           (sub-set1 (mapcar (lambda (prod) (remove term prod :test #'amp-type-iden?))
                             set1))
           (singles (mapcar #'car (remove-if-not #'single-term? sub-set1)))
           (multi-s (remove-if #'single-term? sub-set1)))
      (values singles multi-s set-rest))))

(defun regroup-amp-iter (poly)
  (cond ((null poly) (list '()))
        ((last-one? poly) (list poly))
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
                   (mapcar (lambda (x) (cons '* x)) poly)
                   (mapcan (lambda (common-term) (expand (cadr common-term)))
                           commons)))))))








(defun count-fp (no nv))
