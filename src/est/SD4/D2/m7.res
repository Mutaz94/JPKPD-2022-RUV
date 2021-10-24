Sun Oct 24 04:29:48 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1548.92786457116        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1311E+02  4.9623E+01 -2.0872E+01  1.2562E+02  5.8246E+01 -5.0429E+01 -2.3998E+01 -3.0854E+00 -2.3306E+01 -1.1690E+01
            -1.5190E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1562.35763530368        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0357E+00  9.9623E-01  1.0298E+00  9.5585E-01  9.2832E-01  1.4726E+00  1.1019E+00  1.0137E+00  1.1214E+00  1.0371E+00
             1.0468E+00
 PARAMETER:  1.3506E-01  9.6219E-02  1.2940E-01  5.4850E-02  2.5624E-02  4.8703E-01  1.9706E-01  1.1358E-01  2.1461E-01  1.3641E-01
             1.4575E-01
 GRADIENT:   6.2344E+01  2.7063E+00  3.3070E+01 -3.8809E+01 -7.4221E+01  4.8412E+01 -1.0948E+01 -2.9488E+00 -3.1230E+00  4.6455E+00
             8.9991E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1564.93983381973        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      334
 NPARAMETR:  1.0323E+00  1.0032E+00  9.7513E-01  9.4876E-01  9.1445E-01  1.4472E+00  1.4628E+00  1.0864E+00  1.0862E+00  9.6909E-01
             1.0076E+00
 PARAMETER:  1.3179E-01  1.0317E-01  7.4820E-02  4.7398E-02  1.0566E-02  4.6964E-01  4.8038E-01  1.8284E-01  1.8268E-01  6.8605E-02
             1.0760E-01
 GRADIENT:   6.1025E+01  5.8615E+00  2.8103E+01 -4.0995E+01 -7.0401E+01  4.2929E+01  1.0597E+01  1.6083E+00  8.3329E+00  6.0073E+00
            -2.5504E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1571.77904160447        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  9.7863E-01  9.5174E-01  9.9275E-01  1.0044E+00  9.4935E-01  1.2710E+00  1.4325E+00  9.8503E-01  1.0165E+00  9.7532E-01
             1.0208E+00
 PARAMETER:  7.8402E-02  5.0532E-02  9.2728E-02  1.0438E-01  4.8026E-02  3.3981E-01  4.5943E-01  8.4919E-02  1.1634E-01  7.5006E-02
             1.2058E-01
 GRADIENT:   3.6741E+00  6.5355E-01 -6.6798E-01  3.2008E+00  1.8776E+00 -1.8488E+00  7.3169E-02 -3.7904E-02  2.2606E-01 -6.5304E-01
             9.3543E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1571.84537967505        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  9.7703E-01  1.0578E+00  9.0124E-01  9.3213E-01  9.5508E-01  1.2781E+00  1.3408E+00  9.0220E-01  1.0511E+00  9.7198E-01
             1.0205E+00
 PARAMETER:  7.6766E-02  1.5620E-01 -3.9875E-03  2.9722E-02  5.4039E-02  3.4536E-01  3.9326E-01 -2.9233E-03  1.4986E-01  7.1584E-02
             1.2033E-01
 GRADIENT:   2.4467E-02  1.2468E-01  8.9147E-02  2.5498E-01 -3.3649E-01  1.9487E-01 -2.4135E-01 -3.4874E-03 -2.1877E-01  7.7913E-02
             9.3659E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1571.85245429972        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      865             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7933E-01  1.0750E+00  8.8676E-01  9.2018E-01  9.5695E-01  1.2970E+00  1.3307E+00  8.9003E-01  1.0596E+00  9.7068E-01
             1.0204E+00
 PARAMETER:  7.9111E-02  1.7231E-01 -2.0177E-02  1.6817E-02  5.5995E-02  3.6005E-01  3.8569E-01 -1.6504E-02  1.5791E-01  7.0240E-02
             1.2015E-01
 GRADIENT:   4.4529E+02  7.8931E+01  1.5822E+00  7.4380E+01  7.8596E+00  2.9500E+02  2.9325E+01  1.6768E-01  1.3768E+01  7.8467E-01
             1.0776E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1571.86140214311        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1023
 NPARAMETR:  9.7848E-01  1.0751E+00  8.8671E-01  9.2017E-01  9.5691E-01  1.2886E+00  1.3302E+00  8.8997E-01  1.0591E+00  9.7045E-01
             1.0203E+00
 PARAMETER:  7.8243E-02  1.7238E-01 -2.0237E-02  1.6803E-02  5.5955E-02  3.5353E-01  3.8535E-01 -1.6569E-02  1.5746E-01  7.0006E-02
             1.2005E-01
 GRADIENT:   1.8250E+00 -3.4836E-01 -1.4753E-01 -2.1750E-01  2.1842E-01  3.5398E+00  1.2863E-01  3.1918E-02  7.6925E-02  4.2972E-02
             4.1539E-02

0ITERATION NO.:   33    OBJECTIVE VALUE:  -1571.86145608824        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     1120
 NPARAMETR:  9.7848E-01  1.0752E+00  8.8688E-01  9.2019E-01  9.5683E-01  1.2886E+00  1.3303E+00  8.8868E-01  1.0590E+00  9.7030E-01
             1.0202E+00
 PARAMETER:  7.8244E-02  1.7255E-01 -2.0050E-02  1.6828E-02  5.5868E-02  3.5353E-01  3.8544E-01 -1.8022E-02  1.5736E-01  6.9848E-02
             1.2000E-01
 GRADIENT:   1.4894E-03  1.0164E-02  8.5196E-02 -3.4609E-02 -4.8161E-02  2.6977E-04 -1.7340E-02 -8.7767E-03 -6.5508E-03 -8.4839E-03
            -1.2698E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1120
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7515E-04 -6.6837E-03 -3.0108E-02  1.3738E-03 -2.9163E-02
 SE:             2.9903E-02  2.2038E-02  1.1847E-02  2.3265E-02  2.1450E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9533E-01  7.6168E-01  1.1043E-02  9.5291E-01  1.7395E-01

 ETASHRINKSD(%)  1.0000E-10  2.6169E+01  6.0310E+01  2.2060E+01  2.8141E+01
 ETASHRINKVR(%)  1.0000E-10  4.5490E+01  8.4247E+01  3.9254E+01  4.8363E+01
 EBVSHRINKSD(%)  2.7533E-01  2.6088E+01  6.4300E+01  2.2889E+01  2.5726E+01
 EBVSHRINKVR(%)  5.4990E-01  4.5370E+01  8.7255E+01  4.0539E+01  4.4834E+01
 RELATIVEINF(%)  9.9167E+01  3.2337E+00  1.5457E+00  3.7320E+00  9.2230E+00
 EPSSHRINKSD(%)  4.5336E+01
 EPSSHRINKVR(%)  7.0118E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1571.8614560882390     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -836.71062952450086     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1571.861       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  1.08E+00  8.87E-01  9.20E-01  9.57E-01  1.29E+00  1.33E+00  8.89E-01  1.06E+00  9.70E-01  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       33.244
Stop Time:
Sun Oct 24 04:29:56 CDT 2021
