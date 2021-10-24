Sun Oct 24 01:20:51 CDT 2021
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
$DATA ../../../../data/SD3/D2/dat71.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1929.37138119589        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3219E+02 -4.5539E+01  6.3225E+01 -9.7817E+01  1.3067E+01 -1.8967E+01 -8.4001E+01 -1.9486E+01 -7.3319E+01 -4.2486E+01
            -2.0631E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1970.58638547873        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  9.6307E-01  1.0992E+00  7.0935E-01  9.2828E-01  9.2062E-01  9.8551E-01  1.6986E+00  1.0510E+00  1.1915E+00  1.2045E+00
             9.5444E-01
 PARAMETER:  6.2374E-02  1.9455E-01 -2.4341E-01  2.5573E-02  1.7297E-02  8.5408E-02  6.2980E-01  1.4977E-01  2.7524E-01  2.8606E-01
             5.3366E-02
 GRADIENT:  -3.2539E+01 -4.6663E+01  2.9903E+00 -6.7548E+01 -1.2241E+01 -9.5457E+01  1.3489E-01  1.4777E+01  2.4832E+01  2.7029E+01
            -2.0435E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1992.02311235494        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      289
 NPARAMETR:  9.7493E-01  9.1370E-01  6.3767E-01  1.1016E+00  7.3372E-01  1.2761E+00  2.3335E+00  7.0572E-01  1.0135E+00  8.5698E-01
             9.8790E-01
 PARAMETER:  7.4614E-02  9.7474E-03 -3.4994E-01  1.9679E-01 -2.0962E-01  3.4379E-01  9.4738E-01 -2.4854E-01  1.1341E-01 -5.4335E-02
             8.7824E-02
 GRADIENT:  -1.2405E+00  1.7993E+01  6.4424E+00 -2.2528E+01 -1.9284E+01  2.5697E+01  2.8276E+01  7.7811E+00  1.5825E+01  9.8095E+00
             1.4137E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1997.16957250165        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      467
 NPARAMETR:  9.7531E-01  8.1296E-01  5.9332E-01  1.1520E+00  6.8641E-01  1.1883E+00  2.2603E+00  4.6438E-01  9.8780E-01  8.1512E-01
             9.6837E-01
 PARAMETER:  7.5002E-02 -1.0708E-01 -4.2202E-01  2.4153E-01 -2.7628E-01  2.7255E-01  9.1551E-01 -6.6705E-01  8.7725E-02 -1.0442E-01
             6.7855E-02
 GRADIENT:  -4.4404E-01  3.8539E+00 -1.4233E+01  1.5670E+01  1.1773E+01 -2.1138E+00  1.1359E+01  2.4470E+00  1.1375E+01  4.5890E+00
            -7.9808E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1998.23172844229        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      642
 NPARAMETR:  9.7532E-01  8.3607E-01  5.8046E-01  1.1266E+00  6.8750E-01  1.1930E+00  2.1052E+00  2.7854E-01  9.5871E-01  8.2087E-01
             9.7255E-01
 PARAMETER:  7.5015E-02 -7.9040E-02 -4.4393E-01  2.1917E-01 -2.7469E-01  2.7644E-01  8.4442E-01 -1.1782E+00  5.7835E-02 -9.7385E-02
             7.2162E-02
 GRADIENT:  -1.1405E+00 -4.0460E-01 -1.8990E+00 -1.5124E-02 -2.4912E-02 -6.4337E-01  5.9423E-01  2.9151E-01  9.7364E-01  1.1721E+00
            -4.2601E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1998.28281880405        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  9.7616E-01  8.5251E-01  5.7911E-01  1.1178E+00  6.9438E-01  1.1950E+00  2.0777E+00  1.5684E-01  9.5877E-01  8.3299E-01
             9.7488E-01
 PARAMETER:  7.5869E-02 -5.9566E-02 -4.4626E-01  2.1140E-01 -2.6473E-01  2.7811E-01  8.3125E-01 -1.7525E+00  5.7896E-02 -8.2731E-02
             7.4557E-02
 GRADIENT:  -9.3621E-03 -1.1547E-01 -2.0420E-01 -2.3766E-02  1.4479E-01 -2.7978E-02  3.5935E-02  7.8183E-03  2.0263E-02  6.7744E-02
             1.8994E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1998.28828478975        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      980             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7695E-01  8.5281E-01  5.7920E-01  1.1172E+00  6.9419E-01  1.2003E+00  2.0849E+00  1.4623E-01  9.5852E-01  8.3249E-01
             9.7487E-01
 PARAMETER:  7.6677E-02 -5.9222E-02 -4.4610E-01  2.1083E-01 -2.6502E-01  2.8256E-01  8.3470E-01 -1.8226E+00  5.7631E-02 -8.3336E-02
             7.4550E-02
 GRADIENT:   3.7625E+02  2.6186E+01  2.3801E+01  1.8507E+02  3.4370E+01  1.8799E+02  8.6173E+01  1.5924E-01  7.9068E+00  7.2325E-01
             9.0656E-01

0ITERATION NO.:   33    OBJECTIVE VALUE:  -1998.28861058996        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1072
 NPARAMETR:  9.7709E-01  8.5269E-01  5.7907E-01  1.1176E+00  6.9436E-01  1.1996E+00  2.0842E+00  1.4892E-01  9.5857E-01  8.3295E-01
             9.7491E-01
 PARAMETER:  7.6824E-02 -5.9356E-02 -4.4633E-01  2.1122E-01 -2.6477E-01  2.8200E-01  8.3440E-01 -1.8044E+00  5.7687E-02 -8.2785E-02
             7.4586E-02
 GRADIENT:   2.3752E-01  1.7364E-01  6.8091E-02  3.4722E-01 -1.4217E-02 -2.2171E-01 -1.0174E-01 -1.1608E-03  3.2548E-02  1.8736E-03
            -1.2549E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1072
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6645E-04  1.5085E-02 -9.0429E-03 -1.2405E-02  7.0079E-03
 SE:             2.9945E-02  2.4067E-02  3.3436E-03  2.5509E-02  2.1086E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8491E-01  5.3080E-01  6.8394E-03  6.2676E-01  7.3963E-01

 ETASHRINKSD(%)  1.0000E-10  1.9371E+01  8.8799E+01  1.4541E+01  2.9359E+01
 ETASHRINKVR(%)  1.0000E-10  3.4990E+01  9.8745E+01  2.6968E+01  5.0098E+01
 EBVSHRINKSD(%)  2.3527E-01  1.7225E+01  9.0295E+01  1.4857E+01  3.0171E+01
 EBVSHRINKVR(%)  4.7000E-01  3.1482E+01  9.9058E+01  2.7507E+01  5.1239E+01
 RELATIVEINF(%)  9.9406E+01  1.8110E+01  2.7540E-01  2.0432E+01  8.9460E+00
 EPSSHRINKSD(%)  3.4527E+01
 EPSSHRINKVR(%)  5.7133E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1998.2886105899597     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1079.3500773852870     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1998.289       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  8.53E-01  5.79E-01  1.12E+00  6.94E-01  1.20E+00  2.08E+00  1.49E-01  9.59E-01  8.33E-01  9.75E-01
 


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
 #CPUT: Total CPU Time in Seconds,       35.810
Stop Time:
Sun Oct 24 01:20:59 CDT 2021
