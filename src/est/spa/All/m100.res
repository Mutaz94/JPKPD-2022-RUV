Sat Sep 25 15:17:53 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat100.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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
$OMEGA
0.9 FIX ;     IIV CL
0.9 FIX  ;     IIV V
0.9 FIX ;      IIV KA
$SIGMA  1  FIX;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E+00
 0.0000E+00   0.9000E+00
 0.0000E+00   0.0000E+00   0.9000E+00
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20125.0403938973        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   4.5062E+02 -2.6045E+01  1.6731E+02 -1.7216E+03  1.2785E+03 -1.5604E+02 -4.0875E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -365.081778933073        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  1.3247E+00  1.2846E+00  3.3116E+00  4.1741E-01  3.2224E-01  3.8955E-01  1.5776E+01
 PARAMETER:  3.8117E-01  3.5047E-01  1.2974E+00 -7.7369E-01 -1.0325E+00 -8.4276E-01  2.8585E+00
 GRADIENT:   1.1634E+02  2.0039E+01  4.4186E-01 -3.6058E+01  1.7985E+01  9.6464E-02 -5.9012E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -381.834632960899        NO. OF FUNC. EVALS.:  53
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  1.1007E+00  1.0777E+00  2.0300E+00  3.8046E-01  1.7542E-01  1.9610E-01  1.6021E+01
 PARAMETER:  1.9591E-01  1.7487E-01  8.0801E-01 -8.6638E-01 -1.6406E+00 -1.5291E+00  2.8739E+00
 GRADIENT:  -1.6393E+00 -9.3540E+00 -1.3109E+00 -1.2043E+01  7.3803E+00  4.2002E-01  2.1427E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -386.320311995630        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0457E+00  1.0343E+00  1.7867E+01  3.8565E-01  2.8550E-02  9.8129E-01  1.5704E+01
 PARAMETER:  1.4472E-01  1.3374E-01  2.9830E+00 -8.5282E-01 -3.4561E+00  8.1111E-02  2.8539E+00
 GRADIENT:  -1.4804E+00 -8.9840E+00 -4.0303E-01 -1.6266E+00  8.6422E-02 -7.4434E-03  2.2196E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -386.767332853363        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      212
 NPARAMETR:  1.0480E+00  1.0409E+00  6.8638E+02  3.8889E-01  1.0000E-02  1.8149E+01  1.5719E+01
 PARAMETER:  1.4688E-01  1.4006E-01  6.6314E+00 -8.4445E-01 -5.8983E+00  2.9986E+00  2.8549E+00
 GRADIENT:   4.3938E-01  3.8569E-01 -6.1888E-03 -2.2434E-01  0.0000E+00 -3.6460E-03 -3.6158E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -386.776580063791        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:      281
 NPARAMETR:  1.0485E+00  1.0415E+00  8.6503E+02  3.8897E-01  1.0000E-02  2.1791E+01  1.5740E+01
 PARAMETER:  1.4741E-01  1.4066E-01  6.8628E+00 -8.4426E-01 -6.0577E+00  3.1815E+00  2.8562E+00
 GRADIENT:  -3.5491E-01 -2.4004E-02 -9.0671E-03 -7.1403E-01  0.0000E+00  5.6536E-03 -4.9293E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -386.785307102715        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      361
 NPARAMETR:  1.0490E+00  1.0423E+00  8.6821E+02  3.8915E-01  1.0000E-02  2.1754E+01  1.5771E+01
 PARAMETER:  1.4786E-01  1.4143E-01  6.8664E+00 -8.4378E-01 -6.0577E+00  3.1798E+00  2.8582E+00
 GRADIENT:  -9.1533E-01 -2.1446E-01 -5.1069E-03 -5.4974E-01  0.0000E+00 -2.3439E-03 -3.2200E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -386.799021700259        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      479            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0538E+00  1.0470E+00  1.0755E+03  3.9029E-01  1.0000E-02  2.2439E+01  1.5862E+01
 PARAMETER:  1.5242E-01  1.4595E-01  7.0805E+00 -8.4087E-01 -6.0577E+00  3.2108E+00  2.8639E+00
 GRADIENT:   7.6665E-01  5.4957E-01 -3.9862E-03  6.3099E-01  0.0000E+00 -2.1834E-03  1.9411E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -386.801739233028        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  1.0523E+00  1.0457E+00  1.2447E+04  3.9036E-01  1.0000E-02  3.0861E+01  1.5827E+01
 PARAMETER:  1.5095E-01  1.4473E-01  9.5292E+00 -8.4068E-01 -6.0577E+00  3.5295E+00  2.8617E+00
 GRADIENT:   4.9862E-01  7.4765E-01 -4.9198E-04  6.6592E-01  0.0000E+00 -1.6456E-05  5.8157E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -386.803201690438        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  1.0528E+00  1.0459E+00  1.0302E+04  3.8992E-01  1.0000E-02  2.9883E+01  1.5848E+01
 PARAMETER:  1.5149E-01  1.4489E-01  9.3401E+00 -8.4183E-01 -6.0577E+00  3.4973E+00  2.8630E+00
 GRADIENT:  -1.0257E-01 -7.7202E-02 -5.9837E-04 -1.0473E-01  0.0000E+00 -2.3180E-05 -6.2619E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -386.804243632643        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  1.0536E+00  1.0468E+00  5.0141E+05  3.9009E-01  1.0000E-02  4.3948E+01  1.5865E+01
 PARAMETER:  1.5223E-01  1.4575E-01  1.3225E+01 -8.4139E-01 -6.0577E+00  3.8830E+00  2.8641E+00
 GRADIENT:  -7.4116E-02  9.0942E-02 -1.2696E-05 -1.8560E-02  0.0000E+00 -1.7786E-08 -2.7650E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -386.804260063967        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  1.0537E+00  1.0468E+00  3.3177E+07  3.9010E-01  1.0000E-02  6.8514E+01  1.5866E+01
 PARAMETER:  1.5228E-01  1.4572E-01  1.7417E+01 -8.4135E-01 -6.0577E+00  4.3270E+00  2.8642E+00
 GRADIENT:  -2.0087E-03  2.3394E-03 -1.9214E-07 -4.7420E-04  0.0000E+00  0.0000E+00 -1.1003E-04

0ITERATION NO.:   57    OBJECTIVE VALUE:  -386.804260296514        NO. OF FUNC. EVALS.:  48
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  1.0537E+00  1.0468E+00  9.6109E+07  3.9011E-01  1.0000E-02  7.6432E+01  1.5866E+01
 PARAMETER:  1.5228E-01  1.4572E-01  1.8490E+01 -8.4135E-01 -6.0577E+00  4.4408E+00  2.8642E+00
 GRADIENT:  -1.3006E-03  7.1524E-04  1.9889E-06 -1.1233E-03  0.0000E+00  7.9966E-06  1.9889E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      922
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3694E-02 -4.0613E-03  2.5671E-08
 SE:             8.5878E-02  2.4984E-03  1.9806E-08
 N:                     100         100         100

 P VAL.:         6.1090E-01  1.0405E-01  1.9494E-01

 ETASHRINKSD(%)  9.0207E+00  9.7353E+01  1.0000E+02
 ETASHRINKVR(%)  1.7228E+01  9.9930E+01  1.0000E+02
 EBVSHRINKSD(%)  9.2050E+00  9.7422E+01  1.0000E+02
 EBVSHRINKVR(%)  1.7563E+01  9.9934E+01  1.0000E+02
 RELATIVEINF(%)  9.0210E+00  5.8878E-03  0.0000E+00
 EPSSHRINKSD(%)  2.3084E+00
 EPSSHRINKVR(%)  4.5635E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -386.80426029651403     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       348.34656626722415     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     9.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     1.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -386.804       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.05E+00  1.05E+00  9.70E+07  3.90E-01  1.00E-02  7.68E+01  1.59E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-01
 
 ETA2
+        0.00E+00  9.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.49E-01
 
 ETA2
+        0.00E+00  9.49E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 


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
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        5.44E+02
 
 TH 2
+       -2.73E+02  6.34E+02
 
 TH 3
+       -7.21E-11  1.74E-11  1.05E-20
 
 TH 4
+        4.76E+01 -2.00E+02 -4.61E-12  8.82E+02
 
 TH 5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 6
+        6.02E-05 -1.62E-04  6.35E-15  8.20E-05  0.00E+00  1.99E-07
 
 TH 7
+       -1.00E+01 -1.39E+01  1.73E-14  5.81E+00  0.00E+00  8.58E-07  2.17E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,        9.151
Stop Time:
Sat Sep 25 15:18:09 CDT 2021
