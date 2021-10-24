Sun Oct 24 04:06:03 CDT 2021
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
$DATA ../../../../data/SD4/TD2/dat61.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1700.47462502273        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7227E+02 -2.7527E+01 -1.0471E+01 -9.9062E+00 -1.5380E+01  5.9586E+01  4.2054E+00  1.0996E+01  2.5187E+01  1.9251E+01
            -2.4884E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.34184437874        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0289E+00  1.1285E+00  1.1518E+00  9.8980E-01  1.1275E+00  9.5315E-01  9.8707E-01  9.1644E-01  8.4134E-01  8.8399E-01
             1.1720E+00
 PARAMETER:  1.2845E-01  2.2092E-01  2.4131E-01  8.9752E-02  2.2000E-01  5.2022E-02  8.6990E-02  1.2746E-02 -7.2760E-02 -2.3304E-02
             2.5873E-01
 GRADIENT:   7.6786E-01  8.1378E+00  9.4688E+00  2.2263E-01  1.5869E+00 -9.0633E+00 -5.6150E+00 -1.2785E-01 -1.1028E+01 -1.1365E+01
             2.0679E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1709.54968537350        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0213E+00  9.2312E-01  1.2437E+00  1.1218E+00  1.1021E+00  9.5289E-01  1.0873E+00  3.6644E-01  8.6359E-01  1.0985E+00
             1.1303E+00
 PARAMETER:  1.2109E-01  2.0002E-02  3.1805E-01  2.1495E-01  1.9720E-01  5.1748E-02  1.8372E-01 -9.0393E-01 -4.6657E-02  1.9397E-01
             2.2245E-01
 GRADIENT:  -1.1331E+01  3.4524E+00  1.9735E+00  4.8237E+00  7.1275E-01 -8.0996E+00  1.4218E+00  1.1364E-01  1.4166E+00  2.4817E+00
             1.0899E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1710.23561551433        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0274E+00  9.2969E-01  1.0477E+00  1.1029E+00  1.0053E+00  9.7506E-01  1.1927E+00  3.3744E-01  8.1838E-01  9.6994E-01
             1.0901E+00
 PARAMETER:  1.2704E-01  2.7097E-02  1.4659E-01  1.9794E-01  1.0531E-01  7.4742E-02  2.7619E-01 -9.8636E-01 -1.0043E-01  6.9477E-02
             1.8630E-01
 GRADIENT:   6.3038E-01 -1.1240E+00 -1.5890E-01 -1.6667E+00  4.3506E-01  1.5960E-01  9.6934E-02  1.4592E-02 -3.5078E-01  1.9091E-01
            -6.7042E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1710.28782958424        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0279E+00  1.0260E+00  1.0033E+00  1.0439E+00  1.0243E+00  9.7553E-01  1.1048E+00  3.7383E-01  8.5523E-01  9.5874E-01
             1.0891E+00
 PARAMETER:  1.2755E-01  1.2566E-01  1.0331E-01  1.4292E-01  1.2403E-01  7.5227E-02  1.9962E-01 -8.8396E-01 -5.6382E-02  5.7863E-02
             1.8539E-01
 GRADIENT:  -2.8147E-01  9.1894E-01  2.0476E-01  1.7483E+00 -5.1039E-01 -1.3687E-01  7.9546E-02  3.5512E-02  2.1721E-01  3.4065E-03
            -6.7674E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1710.34570641640        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0293E+00  1.1732E+00  9.2105E-01  9.4894E-01  1.0543E+00  9.7665E-01  9.9756E-01  3.2412E-01  9.1636E-01  9.5293E-01
             1.0898E+00
 PARAMETER:  1.2884E-01  2.5973E-01  1.7760E-02  4.7585E-02  1.5289E-01  7.6370E-02  9.7557E-02 -1.0267E+00  1.2650E-02  5.1790E-02
             1.8604E-01
 GRADIENT:  -3.2840E-02  4.0532E-01 -1.7715E-01  1.3427E+00  1.0910E-01 -2.8680E-01  9.7321E-02  4.7898E-02  4.7195E-01  3.0013E-02
             3.9295E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1710.36279353723        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1078
 NPARAMETR:  1.0300E+00  1.2816E+00  8.5646E-01  8.7944E-01  1.0762E+00  9.7826E-01  9.3415E-01  2.3456E-01  9.6444E-01  9.5288E-01
             1.0896E+00
 PARAMETER:  1.2956E-01  3.4812E-01 -5.4944E-02 -2.8472E-02  1.7340E-01  7.8019E-02  3.1879E-02 -1.3500E+00  6.3794E-02  5.1739E-02
             1.8580E-01
 GRADIENT:  -7.7627E-03  2.5930E+00  9.5562E-01  1.6630E+00 -1.9622E+00  3.4448E-03 -8.7513E-03  1.9204E-02  1.6485E-02  1.1700E-01
            -1.4371E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1710.37252464033        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  1.0313E+00  1.3120E+00  8.3229E-01  8.5752E-01  1.0819E+00  9.7887E-01  9.1911E-01  1.2879E-01  9.7934E-01  9.5108E-01
             1.0896E+00
 PARAMETER:  1.3087E-01  3.7153E-01 -8.3572E-02 -5.3706E-02  1.7874E-01  7.8647E-02  1.5654E-02 -1.9496E+00  7.9121E-02  4.9841E-02
             1.8579E-01
 GRADIENT:   2.6110E+00 -4.7827E-01  3.8940E-01 -4.8683E-01 -1.6928E-01  1.5505E-01  2.1057E-03  4.0479E-03 -4.0038E-02 -1.1873E-01
            -2.0969E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1710.37348473892        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1434
 NPARAMETR:  1.0306E+00  1.3153E+00  8.2865E-01  8.5586E-01  1.0819E+00  9.7879E-01  9.1814E-01  7.0931E-02  9.8026E-01  9.5209E-01
             1.0903E+00
 PARAMETER:  1.3017E-01  3.7406E-01 -8.7963E-02 -5.5653E-02  1.7869E-01  7.8562E-02  1.4598E-02 -2.5460E+00  8.0066E-02  5.0902E-02
             1.8648E-01
 GRADIENT:   8.7900E-01  1.5501E-01  1.0695E-01  4.2250E-01 -1.9122E-02  1.1084E-01  4.7301E-02  1.4135E-03 -3.3457E-02  5.7709E-02
             8.7710E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1710.37476374442        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1615             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0313E+00  1.3150E+00  8.2800E-01  8.5557E-01  1.0816E+00  9.7883E-01  9.1780E-01  1.0000E-02  9.8084E-01  9.5197E-01
             1.0903E+00
 PARAMETER:  1.3081E-01  3.7382E-01 -8.8741E-02 -5.5982E-02  1.7849E-01  7.8600E-02  1.4219E-02 -5.3590E+00  8.0649E-02  5.0783E-02
             1.8642E-01
 GRADIENT:   4.8065E+02  1.7752E+02  1.0859E+00  4.2599E+01  9.7076E+00  4.1384E+01  3.2544E+00  0.0000E+00  3.3389E+00  5.0630E-01
             1.5134E+00

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1710.37476374442        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1670
 NPARAMETR:  1.0313E+00  1.3150E+00  8.2800E-01  8.5557E-01  1.0816E+00  9.7883E-01  9.1780E-01  1.0000E-02  9.8084E-01  9.5197E-01
             1.0903E+00
 PARAMETER:  1.3081E-01  3.7382E-01 -8.8741E-02 -5.5982E-02  1.7849E-01  7.8600E-02  1.4219E-02 -5.3590E+00  8.0649E-02  5.0783E-02
             1.8642E-01
 GRADIENT:   2.3846E+00 -5.6782E-01  4.3427E-02 -9.7385E-02  1.8190E-01  1.2740E-01  1.3856E-02  0.0000E+00  2.4777E-02  2.5651E-02
             2.2004E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1670
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.2818E-04 -1.3881E-02 -3.1461E-04  5.5136E-03 -2.7177E-02
 SE:             2.9811E-02  2.1555E-02  1.3737E-04  2.3036E-02  2.2480E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9122E-01  5.1958E-01  2.2011E-02  8.1084E-01  2.2668E-01

 ETASHRINKSD(%)  1.2916E-01  2.7789E+01  9.9540E+01  2.2827E+01  2.4690E+01
 ETASHRINKVR(%)  2.5816E-01  4.7856E+01  9.9998E+01  4.0443E+01  4.3285E+01
 EBVSHRINKSD(%)  5.1201E-01  2.7105E+01  9.9577E+01  2.3855E+01  2.3335E+01
 EBVSHRINKVR(%)  1.0214E+00  4.6864E+01  9.9998E+01  4.2019E+01  4.1225E+01
 RELATIVEINF(%)  9.8544E+01  1.8557E+00  1.5507E-04  2.2241E+00  7.3265E+00
 EPSSHRINKSD(%)  4.1851E+01
 EPSSHRINKVR(%)  6.6187E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1710.3747637444230     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -975.22393718068486     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1710.375       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.31E+00  8.28E-01  8.56E-01  1.08E+00  9.79E-01  9.18E-01  1.00E-02  9.81E-01  9.52E-01  1.09E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       48.686
Stop Time:
Sun Oct 24 04:06:14 CDT 2021
