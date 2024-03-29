$PROB template control stream 
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained 
; 	   	variability in nonlinear mixed-effect approach 
; Model: 	One-compartment model with linear elimination 
; Estim:	First-order conditional est. with interaction 
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/28/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID 
$DATA ../../data/SD1/B/dat1.csv ignore=@ 
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2 
$PK

	ET1 = EXP(ETA(1)*THETA(4))
	ET2 = EXP(ETA(2)*THETA(5))
	ET3 = EXP(ETA(3)*THETA(6))
	
	
	CL = 5.0 * THETA(1) * ET1
	V = 85  * THETA(2) * ET2 
	KA = 0.7 * THETA(3) * ET3
	
	SC = V
$ERROR
	CVERR 	= 0.05
	W  	= THETA(7)*F*CVERR
	Y  	= F + W * ERR(1)
$THETA  
	 (0,1) ; tvCL
	 (0,1) ; tvV
	 (0,1) ; tvKA
	 (0,1) ; tvCL
	 (0,1) ; tvV
	 (0,1) ; tvK
	 (0,1) ; RUV
$OMEGA (0.09 FIX)x3  
$SIGMA  1  FIX;        [P]
