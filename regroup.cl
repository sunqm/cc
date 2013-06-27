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

;(defun find-most-n-common-term* (n lst-term))

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

;;; given a list of products to be sumed up, regroup them to be
;;; a list of polynomials; each polynomial is sth like (+ (* a b (+ c ..)))
(defun regroup-iter (poly)
  (if (null poly) '()
      (let* ((match (find-most-n-common-term 1 poly))
             (count (caar match))
             (obj (cadar match)))
        (if (eql count 1)
            poly
            (multiple-value-bind (singles multi-s set-rest)
                (subgroup-terms obj poly)
              `((* ,obj (+ ,@singles ,@(regroup multi-s)))
                ,@(regroup set-rest)))))))
(defun regroup (poly)
  (let ((res (regroup-iter poly)))
    (if (last-one? res)
        (car res)
        (cons '+ (regroup-iter poly)))))












(defun count-fp (no nv))
