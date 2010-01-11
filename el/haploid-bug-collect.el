;;; haploid-bug-collect.el --- Collect and file bug for haploid library

;; Copyright (C) 2010  Joel James Adamson

;; Author: Joel James Adamson <adamsonj@email.unc.edu>
;; Keywords: mail, tools

;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:

;; Scan a mail buffer and collect information from a bug report;
;; insert this into the chosen working copy's BUGS file (in org-mode)
;;
;; This library provides two variables and one command
;;
;; To collect bug data and insert it into your working copy, use
;;
;; M-x haploid-bug-collect

;;; Code:
(defgroup haploid nil
  "Haploid genetics library project management."
  :group 'programming
  :prefix 'haploid-bug-collect)

(defcustom haploid-bug-collect-working-bugs-file nil
  "Your working copy's BUGS file (full path name)"
  :group 'haploid
  :type '(file :must-match t))

(defcustom haploid-bug-collect-hook nil
  "Hook to call after collecting bug information and inserting
bug report"
  :group 'haploid
  :type 'hook)

(provide 'haploid-bug-collect)
;;; haploid-bug-collect.el ends here
