;;; -*- Mode: Common Lisp; -*-
;;;
;;; 

(load "utility.cl")

; regroup polynomial
(defun sum? (poly)
  (eql '+ (car poly)))
(defun prod? (poly)
  (eql '* (car poly)))
(defun remove-prod-operator (term)
  (if (eql '* (car term))
      (cdr term)
      term))
(defun intersectionp (&rest lsts)
  (let ((lst1 (car lsts))
        (lst-rest (cdr lsts)))
    (some (lambda (item)
            (every (lambda (lst) (member item lst :test #'equal))
                   lst-rest))
          lst1)))

;;; find unique terms which are reasonable for regrouping.
;;; define a filter function which should consider memory usage
;;; (restriction on the number of open lines)
(defun find-uniq-term (lst-term &key (filter #'remove-prod-operator))
  (remove-duplicates (mapcan filter lst-term) ; first flatten the list
                     :test #'equal))

(defun find-most-n-common-term (n lst-term)
  (labels ((count-obj (obj)
             (count-if (lambda (term) (member obj term :test #'equal))
                       lst-term))
           (keep-terms-with (keys) ;find the exclusive terms during previous search
             (lambda (term)
               (if (subsetp keys term)
                   (remove-if (lambda (x) (member x keys :test #'equal))
                              (remove-prod-operator term)))))
           (find-most (n)
             (if (eql n 1)
                 (multiple-value-bind (obj count)
                     (most #'count-obj (find-uniq-term lst-term))
                   (list (list count obj)))
                 (let* ((matches (find-most (1- n)))
                        (keys (mapcar #'cadr matches))
                        (objs (find-uniq-term lst-term :filter (keep-terms-with keys))))
                   (if objs
                       (multiple-value-bind (obj count) (most #'count-obj objs)
                         (cons (list count obj) matches))
                       matches)))))
    (reverse (find-most n))))

;;; return two sets: the first set includes obj and the second set collects the rest
(defun subgroup-terms (obj lst-term)
  (flet ((contain-obj (term)
           (member obj term :test #'equal))
         (single-term? (term)
           (last-one? (remove-prod-operator term))))
    (let* ((set1 (remove-if-not #'contain-obj lst-term))
           (set-rest (remove-if #'contain-obj lst-term))
           (sub-set1 (mapcar (lambda (term) (remove obj term :test #'equal))
                             set1))
           (singles (mapcar #'cadr (remove-if-not #'single-term? sub-set1)))
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
        (t (flet ((expand (common-term)
                    (let ((count (car common-term))
                          (obj (cadr common-term)))
                      (if (eql count 1)
                          '() ; don't regroup since no common terms are found
                          (multiple-value-bind (singles multi-s set-rest)
                              (subgroup-terms obj poly)
                            (flet ((comb (regp-mults regp-rest)
                                     `((* ,obj (+ ,@singles ,@regp-mults))
                                       ,@regp-rest)))
                              (mapnest #'comb
                                       (regroup-iter multi-s)
                                       (regroup-iter set-rest))))))))
             (mapcan #'expand (find-most-n-common-term *common-term-search-trial* poly))))))
(defun regroup (poly)
  (mapcar (lambda (res) (cons '+ res))
          (regroup-iter poly)))












(defun count-fp (no nv))
