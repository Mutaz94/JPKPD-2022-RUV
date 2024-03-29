$PROB template control stream 
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained 
; 	   	variability components in nonlinear mixed-effect approach 
; Model: 	Two-compartment model with linear elimination 
; Estim:		First-order conditional est. with interaction 
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID 
$DATA <DATA> ignore=@ 
$SUBR ADVAN4 TRANS4
$EST MET=1 NOA MAX=10000 PRINT=5 INT NSIG=2 
$PK

;; 	Take SQRT since the values are in variance term

	ET1 = EXP(ETA(1)*THETA(6))
	ET2 = EXP(ETA(2)*THETA(7))
	ET3 = EXP(ETA(3)*THETA(8))
	ET4 = EXP(ETA(4)*THETA(9))
	ET5 = EXP(ETA(5)*THETA(10))
	
;;	True values 

	TCL = 5.0
	TV2 = 35.0
	TQ  = 50
	TV3 = 50 
	TKA = 0.7 
	CVERR = 0.05



	CL = TCL  * THETA(1) * ET1
	V2 = TV2  * THETA(2) * ET2 
	Q  = TQ   * THETA(3) * ET3
	V3 = TV3  * THETA(4) * ET4
	KA = TKA  * THETA(5) * ET5	
	SC = V2 

$ERROR
	
	W = THETA(11)*F*CVERR 
	Y 	= F + W*ERR(1)

$THETA  
	 (0,1) ; CL
	 (0,1) ; V2
	 (0,1) ; Q
	 (0,1) ; V3
	 (0,1) ; KA
	 (0,1) ; IIVCL
	 (0,1) ; IIVV2
	 (0,1) ; IIVQ
	 (0,1) ; IIVV3
	 (0,1) ; IIVKA
	 (0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5   
$SIGMA  1 FIX ;        
