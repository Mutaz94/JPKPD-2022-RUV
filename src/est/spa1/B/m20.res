Wed Sep 29 20:56:09 CDT 2021
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
$DATA ../../../../data/spa1/B/dat20.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1546.56002514075        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8844E+02  1.5049E+00  3.2352E+01  4.2195E+01  8.3779E+01  3.2955E+01  1.4661E+00 -2.1080E+02 -3.7023E+01 -1.0185E+01
            -8.1738E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2059.51023168269        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:      144             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5670E-01  1.0127E+00  1.0327E+00  9.9134E-01  9.5457E-01  1.0683E+00  9.9221E-01  1.1519E+00  1.0047E+00  8.6965E-01
             1.4191E+00
 PARAMETER:  5.5735E-02  1.1264E-01  1.3213E-01  9.1305E-02  5.3503E-02  1.6611E-01  9.2178E-02  2.4138E-01  1.0469E-01 -3.9662E-02
             4.5002E-01
 GRADIENT:   1.9503E+02  1.6100E+01 -4.5719E+00  3.1866E+01 -1.5522E+00  9.1697E+01  1.3539E-01  1.0934E+01  7.0971E+00  2.2785E+00
             2.0609E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2060.79905816372        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.7613E-01  1.0127E+00  1.0845E+00  1.0028E+00  9.7420E-01  1.0119E+00  1.0781E+00  1.1519E+00  9.7648E-01  8.6104E-01
             1.4191E+00
 PARAMETER:  7.5844E-02  1.1265E-01  1.8109E-01  1.0279E-01  7.3859E-02  1.1183E-01  1.7518E-01  2.4137E-01  7.6194E-02 -4.9616E-02
             4.5001E-01
 GRADIENT:   3.7508E+00  1.6655E+00 -4.6298E-01 -9.2845E-01 -3.8351E+00  9.0421E-02  9.3334E-01  8.3805E+00  4.7572E-01 -1.4976E+00
             2.0024E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2060.99363694950        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      477
 NPARAMETR:  9.7432E-01  1.0127E+00  1.1714E+00  1.0061E+00  1.0154E+00  1.0123E+00  1.0146E+00  1.1519E+00  9.9220E-01  9.4079E-01
             1.4190E+00
 PARAMETER:  7.3981E-02  1.1264E-01  2.5822E-01  1.0611E-01  1.1532E-01  1.1223E-01  1.1451E-01  2.4137E-01  9.2166E-02  3.8965E-02
             4.4993E-01
 GRADIENT:   5.2267E-01 -1.8141E-01 -6.7014E-02 -1.7171E+00  6.2657E-01  2.2317E-01 -2.8806E-01  6.0132E+00 -4.4636E-01  1.5208E-02
             2.0165E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2062.39320533156        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      666             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7424E-01  9.8770E-01  1.1830E+00  1.0067E+00  1.0183E+00  1.0131E+00  9.1160E-01  1.1518E+00  1.0193E+00  9.5034E-01
             1.4067E+00
 PARAMETER:  7.3900E-02  8.7627E-02  2.6807E-01  1.0663E-01  1.1811E-01  1.1304E-01  7.4513E-03  2.4129E-01  1.1916E-01  4.9060E-02
             4.4122E-01
 GRADIENT:   2.3733E+02  5.7606E+00 -2.5094E+00  2.5831E+01  1.5079E+01  4.8472E+01 -1.4012E-01  6.2021E+00  6.2654E+00 -1.0734E+00
             2.0168E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2062.74992088604        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      842
 NPARAMETR:  9.7415E-01  1.0015E+00  1.1830E+00  1.0150E+00  1.0145E+00  1.0123E+00  1.0302E+00  1.1518E+00  9.8736E-01  9.4613E-01
             1.4066E+00
 PARAMETER:  7.3809E-02  1.0152E-01  2.6808E-01  1.1492E-01  1.1444E-01  1.1222E-01  1.2977E-01  2.4129E-01  8.7280E-02  4.4628E-02
             4.4120E-01
 GRADIENT:   6.6369E-01  1.7405E+00  5.1819E-01  4.5648E-01 -4.6925E-01  2.1573E-01  1.0016E-01  5.8836E+00  3.0961E-02  1.8465E-01
             1.9801E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2062.92795733324        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  9.7281E-01  8.3422E-01  1.1833E+00  1.1193E+00  9.4456E-01  1.0113E+00  1.2291E+00  1.1518E+00  8.9746E-01  8.6873E-01
             1.4052E+00
 PARAMETER:  7.2433E-02 -8.1264E-02  2.6829E-01  2.1274E-01  4.2960E-02  1.1123E-01  3.0627E-01  2.4128E-01 -8.1846E-03 -4.0725E-02
             4.4018E-01
 GRADIENT:   4.8881E-01  1.5786E+00 -3.3812E+00  3.5063E+00 -1.6618E+00  3.8541E-01 -3.8932E-01  8.7366E+00 -1.1852E+00 -1.0238E-01
             1.9682E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2063.00907209453        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  9.7083E-01  6.5493E-01  1.1838E+00  1.2249E+00  8.8042E-01  1.0107E+00  1.5354E+00  1.1517E+00  8.3144E-01  8.0105E-01
             1.4021E+00
 PARAMETER:  7.0394E-02 -3.2323E-01  2.6875E-01  3.0286E-01 -2.7356E-02  1.1066E-01  5.2882E-01  2.4127E-01 -8.4600E-02 -1.2184E-01
             4.3798E-01
 GRADIENT:   3.1400E-01 -4.3396E-01 -7.8166E+00  1.4166E+00 -1.2309E+00  9.0359E-01  4.2718E-02  1.1191E+01 -2.3076E-01 -3.0834E-02
             1.9479E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2063.15343871566        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1378
 NPARAMETR:  9.6740E-01  4.4141E-01  1.1851E+00  1.3477E+00  8.1385E-01  1.0117E+00  2.0917E+00  1.1517E+00  7.6183E-01  7.3575E-01
             1.3952E+00
 PARAMETER:  6.6853E-02 -7.1778E-01  2.6979E-01  3.9838E-01 -1.0598E-01  1.1163E-01  8.3800E-01  2.4123E-01 -1.7203E-01 -2.0686E-01
             4.3301E-01
 GRADIENT:  -5.9353E-01 -5.6139E-01 -1.0784E+01 -1.2386E-01 -4.4232E+00  2.3048E+00 -4.8131E-01  1.2076E+01 -3.4770E+00 -5.8756E-02
             1.9019E+02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -2063.16585737778        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1444
 NPARAMETR:  9.6732E-01  4.3286E-01  1.1852E+00  1.3525E+00  8.1174E-01  1.0109E+00  2.1209E+00  1.1516E+00  7.6132E-01  7.3377E-01
             1.3949E+00
 PARAMETER:  6.6719E-02 -7.3763E-01  2.6985E-01  4.0207E-01 -1.0887E-01  1.1144E-01  8.5204E-01  2.4123E-01 -1.7444E-01 -2.0961E-01
             4.3272E-01
 GRADIENT:  -9.8468E-01 -3.8696E-01 -3.0964E+04  4.1552E+04 -3.9734E+00  1.8630E+00  1.6853E+04  3.4603E+04 -3.4553E+00 -6.1574E-02
            -3.8648E+04
 NUMSIGDIG:         2.0         2.2         2.3         2.3         1.3         1.1         2.4         2.3         0.8         2.4
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1444
 NO. OF SIG. DIGITS IN FINAL EST.:  0.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6686E-03  2.0134E-02 -2.2785E-02 -1.8878E-02 -2.0211E-02
 SE:             2.9621E-02  1.7372E-02  1.5947E-02  2.5191E-02  1.7879E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5508E-01  2.4645E-01  1.5308E-01  4.5363E-01  2.5828E-01

 ETASHRINKSD(%)  7.6536E-01  4.1802E+01  4.6574E+01  1.5607E+01  4.0105E+01
 ETASHRINKVR(%)  1.5249E+00  6.6130E+01  7.1457E+01  2.8778E+01  6.4126E+01
 EBVSHRINKSD(%)  6.3733E-01  4.6141E+01  4.2666E+01  1.5125E+01  3.8502E+01
 EBVSHRINKVR(%)  1.2706E+00  7.0992E+01  6.7128E+01  2.7962E+01  6.2179E+01
 RELATIVEINF(%)  9.7851E+01  2.8692E+00  5.5659E+00  8.1553E+00  6.6394E+00
 EPSSHRINKSD(%)  4.7853E+01
 EPSSHRINKVR(%)  7.2807E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2063.1658573777777     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1144.2273241731050     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2063.166       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.67E-01  4.33E-01  1.19E+00  1.35E+00  8.11E-01  1.01E+00  2.12E+00  1.15E+00  7.60E-01  7.34E-01  1.39E+00
 


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
+        1.15E+03
 
 TH 2
+       -2.11E+01  4.68E+06
 
 TH 3
+       -9.54E+01 -4.09E+06  4.08E+06
 
 TH 4
+        5.19E+01 -3.94E+02  3.39E+05  1.41E+06
 
 TH 5
+        4.43E+00 -2.95E+02 -2.08E+06  8.32E+02  1.20E+03
 
 TH 6
+        2.53E+00 -3.87E+00 -1.74E+02  1.00E+02 -9.56E-01  1.87E+02
 
 TH 7
+        3.04E+01  7.24E+05 -7.22E+05 -3.90E+01 -3.68E+05  5.57E+01  1.28E+05
 
 TH 8
+        1.09E+02 -1.60E+03  4.63E+03 -3.92E+05  1.53E+03  2.00E+02 -7.85E+01  5.40E+06
 
 TH 9
+        3.39E+00  1.13E+07  1.08E+03 -6.50E+02 -3.57E+07 -2.80E+07  1.74E+06 -1.23E+03  2.71E+07
 
 TH10
+        1.48E+00  1.24E+01 -8.49E+06 -4.81E+02 -8.72E+01  2.92E-01 -1.36E+02 -8.96E+02 -3.40E+00  6.47E+01
 
 TH11
+       -6.06E+01  7.31E+02 -2.16E+03  1.34E+03 -7.10E+02 -9.09E+01  3.77E+01  7.31E+03  5.83E+02  4.40E+02  1.16E+06
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       30.898
Stop Time:
Wed Sep 29 20:56:41 CDT 2021
