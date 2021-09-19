Sat Sep 18 13:16:14 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/S2/dat14.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m14.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1706.34787744647        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5112E+01 -1.0177E+02 -6.1220E+01 -5.9159E+01  7.4374E+01 -3.3169E+00  1.6780E+00  1.3985E+01  3.0488E+01  9.9216E+00
             1.0238E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1715.97943288440        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9129E-01  1.1342E+00  1.1473E+00  1.0032E+00  1.0518E+00  1.0099E+00  9.5451E-01  8.9152E-01  8.0905E-01  9.3043E-01
             1.0036E+00
 PARAMETER:  9.1248E-02  2.2595E-01  2.3739E-01  1.0321E-01  1.5052E-01  1.0987E-01  5.3448E-02 -1.4826E-02 -1.1189E-01  2.7890E-02
             1.0355E-01
 GRADIENT:  -4.4366E+00  3.7730E+01  3.1782E+00  4.3309E+01  1.6023E+01  6.5312E-01 -5.2414E+00  7.7821E-01 -1.2703E+01 -1.3900E+01
            -2.2336E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1718.30102113674        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0052E+00  1.0445E+00  1.0003E+00  1.0449E+00  9.8015E-01  1.0425E+00  9.9800E-01  2.9347E-01  9.0391E-01  9.8887E-01
             1.0106E+00
 PARAMETER:  1.0515E-01  1.4351E-01  1.0030E-01  1.4390E-01  7.9953E-02  1.4164E-01  9.8003E-02 -1.1260E+00 -1.0290E-03  8.8807E-02
             1.1056E-01
 GRADIENT:   2.7237E+01  4.1411E-01 -2.0328E+01  4.2966E+01  3.6906E+01  1.4820E+01  3.7704E-01  4.1634E-01  1.2391E+01 -1.8937E+00
             1.0441E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1719.19864191407        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9679E-01  1.0592E+00  9.6304E-01  1.0189E+00  9.5520E-01  1.0134E+00  1.0316E+00  3.4137E-01  8.5270E-01  9.6270E-01
             9.9252E-01
 PARAMETER:  9.6787E-02  1.5748E-01  6.2340E-02  1.1870E-01  5.4165E-02  1.1330E-01  1.3115E-01 -9.7479E-01 -5.9353E-02  6.1988E-02
             9.2496E-02
 GRADIENT:   7.4108E+00 -1.2411E+00 -7.2578E+00  1.0417E+01  1.2580E+01  2.3595E+00  9.3874E-01  6.6084E-01  2.9391E+00  3.0013E-01
             2.7111E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1719.21133324289        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9449E-01  1.0600E+00  9.4073E-01  1.0143E+00  9.4060E-01  1.0102E+00  1.0330E+00  2.8545E-01  8.4476E-01  9.5141E-01
             9.8891E-01
 PARAMETER:  9.4476E-02  1.5823E-01  3.8897E-02  1.1415E-01  3.8758E-02  1.1019E-01  1.3251E-01 -1.1537E+00 -6.8708E-02  5.0192E-02
             8.8850E-02
 GRADIENT:   2.1850E+00 -1.3540E+00 -3.3007E+00  3.4790E+00  4.6774E+00  7.9672E-01  3.1116E-01  4.6639E-01  1.3788E+00  5.1147E-01
             1.2437E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1719.22290433344        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.9333E-01  1.0617E+00  9.2169E-01  1.0106E+00  9.2985E-01  1.0086E+00  1.0349E+00  2.1080E-01  8.3926E-01  9.4285E-01
             9.8655E-01
 PARAMETER:  9.3309E-02  1.5983E-01  1.8457E-02  1.1055E-01  2.7264E-02  1.0852E-01  1.3434E-01 -1.4568E+00 -7.5237E-02  4.1154E-02
             8.6463E-02
 GRADIENT:  -5.9741E-01 -9.8355E-01 -7.8847E-01 -3.6335E-01  3.7711E-02 -7.8514E-02 -4.8049E-02  2.5281E-01  2.9603E-01  4.9342E-01
             2.7647E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1719.33592492624        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.9356E-01  1.0644E+00  9.0924E-01  1.0083E+00  9.2501E-01  1.0086E+00  1.0367E+00  1.6612E-02  8.3780E-01  9.3856E-01
             9.8595E-01
 PARAMETER:  9.3539E-02  1.6241E-01  4.8562E-03  1.0829E-01  2.2051E-02  1.0858E-01  1.3606E-01 -3.9976E+00 -7.6970E-02  3.6595E-02
             8.5855E-02
 GRADIENT:  -2.8349E-01  2.7487E-03  2.8290E-01 -3.2389E-01 -4.6854E-01 -7.8643E-02 -6.9656E-02  1.5276E-03 -8.7772E-02 -5.2468E-02
            -7.0426E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1719.91300735401        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  1.0165E+00  1.1474E+00  8.9829E-01  9.6499E-01  9.5650E-01  1.0238E+00  9.7689E-01  1.0000E-02  8.7407E-01  9.5385E-01
             9.8914E-01
 PARAMETER:  1.1634E-01  2.3751E-01 -7.2667E-03  6.4358E-02  5.5525E-02  1.2348E-01  7.6615E-02 -9.5248E+00 -3.4599E-02  5.2755E-02
             8.9079E-02
 GRADIENT:   4.2716E+00  3.1110E+00  8.4750E-01  4.1572E+00 -6.7605E-01  1.0256E+00  1.1811E-01  0.0000E+00 -1.9618E-01 -5.1576E-01
             2.3921E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1719.97251919579        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  1.0151E+00  1.2639E+00  8.4846E-01  8.8820E-01  9.8902E-01  1.0218E+00  8.9874E-01  1.0000E-02  9.3234E-01  9.6258E-01
             9.8870E-01
 PARAMETER:  1.1498E-01  3.3422E-01 -6.4337E-02 -1.8564E-02  8.8962E-02  1.2156E-01 -6.7614E-03 -9.9235E+00  2.9942E-02  6.1863E-02
             8.8640E-02
 GRADIENT:   3.2265E-01  2.0267E-01  9.8483E-02  8.5440E-02 -2.2032E-01  2.1355E-02 -2.5388E-02  0.0000E+00 -8.7113E-02  4.4573E-02
            -4.2829E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1719.97272771904        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  1.0150E+00  1.2712E+00  8.4528E-01  8.8344E-01  9.9120E-01  1.0218E+00  8.9454E-01  1.0000E-02  9.3672E-01  9.6274E-01
             9.8885E-01
 PARAMETER:  1.1485E-01  3.3994E-01 -6.8091E-02 -2.3934E-02  9.1164E-02  1.2155E-01 -1.1450E-02 -9.9944E+00  3.4634E-02  6.2033E-02
             8.8788E-02
 GRADIENT:   1.6468E-03 -3.5426E-03  7.3417E-03 -1.0362E-02 -8.0627E-03  1.9777E-04  2.5768E-05  0.0000E+00  2.2690E-03  3.1037E-04
             1.1367E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      908
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.0820E-04 -1.4338E-02 -3.2873E-04  5.8874E-03 -2.3034E-02
 SE:             2.9838E-02  2.1051E-02  1.3813E-04  2.3350E-02  2.3767E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9443E-01  4.9581E-01  1.7321E-02  8.0094E-01  3.3247E-01

 ETASHRINKSD(%)  3.9040E-02  2.9478E+01  9.9537E+01  2.1774E+01  2.0377E+01
 ETASHRINKVR(%)  7.8064E-02  5.0266E+01  9.9998E+01  3.8807E+01  3.6602E+01
 EBVSHRINKSD(%)  3.9738E-01  2.9178E+01  9.9568E+01  2.2635E+01  1.8631E+01
 EBVSHRINKVR(%)  7.9318E-01  4.9842E+01  9.9998E+01  4.0146E+01  3.3790E+01
 RELATIVEINF(%)  9.8862E+01  1.5850E+00  1.6655E-04  2.0825E+00  7.8062E+00
 EPSSHRINKSD(%)  4.3171E+01
 EPSSHRINKVR(%)  6.7704E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.9727277190398     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.82190115530159     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.973       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.27E+00  8.45E-01  8.83E-01  9.91E-01  1.02E+00  8.95E-01  1.00E-02  9.37E-01  9.63E-01  9.89E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.03E+03
 
 TH 2
+       -6.75E+00  4.50E+02
 
 TH 3
+        1.39E+01  1.47E+02  3.19E+02
 
 TH 4
+       -8.49E+00  4.44E+02 -2.08E+02  9.75E+02
 
 TH 5
+       -1.62E+00 -2.63E+02 -4.20E+02  2.36E+02  8.05E+02
 
 TH 6
+       -2.65E+00 -1.17E+00 -3.39E-01 -5.65E+00 -1.83E+00  1.89E+02
 
 TH 7
+        2.87E+00  2.12E+01  5.01E+00 -1.75E+01 -6.23E+00  3.30E-01  7.27E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.63E+00 -2.21E+01 -2.80E+01  3.33E+01  7.35E+00 -1.41E+00  2.68E+01  0.00E+00  8.99E+01
 
 TH10
+       -3.05E+00 -1.03E+01 -3.97E+01 -1.29E+01 -5.70E+01 -4.32E+00  1.92E+01  0.00E+00  1.02E+01  8.97E+01
 
 TH11
+       -9.26E+00 -2.14E+01 -3.51E+01 -1.98E+00  5.25E+00  3.13E-02  8.19E+00  0.00E+00  1.24E+01  2.03E+01  2.24E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       14.872
Stop Time:
Sat Sep 18 13:16:30 CDT 2021
