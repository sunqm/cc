;;;;
;;;; File: utility.cl
;;;; Author
;;;; Description:

(defun last-one? (lst)
  (null (cdr lst)))
(defun last-two? (lst)
  (null (cddr lst)))

(defmacro caddddr (l)
  `(nth 4 ,l))
(defmacro cadddddr (l)
  `(nth 5 ,l))

; see intersection
;(defun (intersect? set1 set2)
;  (cond
;    ((null set1) '())
;    ((member (car set1) set2) t)
;    (t (intersect? (cdr set1) set2))))
;(defun intersect (set1 set2)
;  (cond
;    ((null set1) '())
;    ((member (car set1) set2)
;     (cons (car set1)
;           (intersect (cdr set1) set2)))
;    (t (intersect (cdr set1) set2))))
;(defun intersect (set1 set2)
;  (remove-if-not (lambda (n)
;                   (member n set1))
;                 set2))

;;; based on On.Lisp
(defun group (source n)
  (if (zerop n) (error "zero length"))
  (labels ((rec (source acc)
             (let ((rest (nthcdr n source)))
               (if (consp rest)
                 (rec rest (cons (subseq source 0 n) acc))
                 (nreverse (cons source acc))))))
    (if source (rec source nil) nil)))

(defmacro while (test &body body)
  `(do ()
     ((not ,test))
     ,@body))

(defun last1 (lst)
  (car (last lst)))

(defun flatten (x)
  (labels ((rec (x acc)
             (cond ((null x) acc)
                   ((atom x) (cons x acc))
                   (t (rec (car x) (rec (cdr x) acc))))))
    (rec x nil)))

; recursive remove-if
(defun prune (test tree)
  (labels ((rec (tree acc)
             (cond ((null tree) (nreverse acc))
                   ((consp (car tree))
                    (rec (cdr tree)
                         (cons (rec (car tree) nil) acc)))
                   (t (rec (cdr tree)
                           (if (funcall test (car tree))
                               acc
                               (cons (car tree) acc)))))))
    (rec tree nil)))
(defun remove-if* (test lst)
  (prune test lst))

(defun split-if (fn lst)
  (let ((acc nil))
    (do ((src lst (cdr src)))
      ((or (null src) (funcall fn (car src)))
       (values (nreverse acc) src))
      (push (car src) acc))))

;;; > (most #'length '((a b) (a b c) (a) (e f g)))
;;; (A B C)
;;; 3
(defun most (fn lst)
  (if (null lst)
      (values nil nil)
      (let* ((wins (car lst))
             (max (funcall fn wins)))
        (dolist (obj (cdr lst))
          (let ((score (funcall fn obj)))
            (when (> score max)
              (setq wins obj
                    max score))))
        (values wins max))))

;;; > (mostn #'length '((a b) (a b c) (a) (e f g)))
;;; ((A B C) (E F G))
;;; 3
(defun mostn (fn lst)
  (if (null lst)
      (values nil nil)
      (let ((result (list (car lst)))
            (max (funcall fn (car lst))))
        (dolist (obj (cdr lst))
          (let ((score (funcall fn obj)))
            (cond ((> score max)
                   (setq max    score
                         result (list obj)))
                  ((= score max)
                   (push obj result)))))
        (values (nreverse result) max))))

(defun best (fn lst)
  (if (null lst)
      nil
      (let ((wins (car lst)))
        (dolist (obj (cdr lst))
          (if (funcall fn obj wins)
              (setq wins obj)))
        wins)))

(defun mapa-b (fn a b &optional (step 1))
  (do ((i a (+ i step))
       (result nil))
      ((> i b) (nreverse result))
    (push (funcall fn i) result)))

(defun map0-n (fn n)
  (mapa-b fn 0 n))

(defun map1-n (fn n)
  (mapa-b fn 1 n))

(defun range (n &optional m)
  (cond ((null m)
         (mapa-b #'1+ -1 (- n 2)))
        (t (mapa-b #'1+ (1- n) (- m 2)))))

(defun map-> (fn start test-fn succ-fn)
  (do ((i start (funcall succ-fn i))
       (result nil))
      ((funcall test-fn i) (nreverse result))
    (push (funcall fn i) result)))

;;; > (mkstr pi " pieces of " 'pi)
;;; "3.141592653589793 pieces of PI"
(defun mkstr (&rest args)
  (with-output-to-string (s)
    (dolist (a args) (princ a s))))

;;; > (symb 'ar "Madi" #\L #\L 0)
;;; |ARMadiLL0| ; return the name of variable
(defun symb (&rest args)
  (values (intern (apply #'mkstr args))))

(defun reread (&rest args)
  (values (read-from-string (apply #'mkstr args))))

;;; > (explode 'bomb)
;;; (B O M B)
(defun explode (sym)
  (map 'list #'(lambda (c)
                 (intern (make-string 1
                                      :initial-element c)))
             (symbol-name sym)))

(defun for-each-comps-apply (comps vs f)
  (mapcar (lambda (comp)
              `(,f ,(mapcar (lambda (v) `(,comp ,v)) vs)))
            comps))


;;; functions added for CC/MP
(defun prepend-for-each (a lst)
  (mapcar (lambda (x) (cons a x)) lst))

(defun append-for-each (a lst)
  (mapcar (lambda (x) (append x `(,a))) lst))

(defun check-for-each (assertor lst)
  (cond ((null lst) t)
        ((last-one? lst) (funcall assertor (car lst)))
        (t (and (funcall assertor (car lst))
                (check-for-each assertor (cdr lst))))))

; 0 based
(defun insert-after-0th (a l)
  (cons (car l)
        (cons a (cdr l))))
; 0 based
;(defun insert-after-nth (n a l))

(defun mapnest (f l1 l2)
  (mapcan (lambda (a1)
            (mapcar (lambda (a2)
                      (funcall f a1 a2))
                    l2))
            l1))
