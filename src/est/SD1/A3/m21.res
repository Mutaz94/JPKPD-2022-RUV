Sat Oct 23 14:11:49 CDT 2021
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
$DATA ../../../../data/SD1/A3/dat21.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1010.53847295655        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7394E+02  3.6482E+02  3.9085E+02  2.2163E+01  3.3874E+02  7.0170E+01 -3.3172E+02 -3.1047E+02 -4.1621E+01 -2.2307E+02
            -4.7589E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2807.49514716897        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0563E+00  8.2834E-01  7.8967E-01  1.1555E+00  7.9357E-01  8.1668E-01  1.0286E+00  9.7173E-01  1.0020E+00  8.8680E-01
             2.6707E+00
 PARAMETER:  1.5474E-01 -8.8337E-02 -1.3615E-01  2.4455E-01 -1.3121E-01 -1.0251E-01  1.2820E-01  7.1322E-02  1.0203E-01 -2.0141E-02
             1.0823E+00
 GRADIENT:   1.3537E+02  4.5143E+01  9.1160E+00  5.5397E+01  2.6229E+01 -2.8670E+01  3.4711E+00  1.4794E+01  2.8019E+00  1.3366E+01
             9.4838E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2879.33060574131        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0810E+00  3.3591E-01  2.4113E-01  1.4995E+00  2.7383E-01  7.3532E-01  1.2272E+00  7.4659E-01  1.2564E+00  6.6561E-01
             2.2533E+00
 PARAMETER:  1.7792E-01 -9.9091E-01 -1.3224E+00  5.0513E-01 -1.1952E+00 -2.0745E-01  3.0471E-01 -1.9223E-01  3.2828E-01 -3.0705E-01
             9.1238E-01
 GRADIENT:   2.8282E+02  2.2522E+02  1.6843E+01  4.8183E+02  2.5541E+02 -9.1508E+01 -1.3826E+01  7.3916E+00 -7.4502E+01  2.2869E+01
             9.5042E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2901.87974318912        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  1.0715E+00  3.1922E-01  2.3149E-01  1.5179E+00  2.6290E-01  7.7803E-01  1.3748E+00  7.9830E-01  1.6564E+00  4.4401E-01
             2.2477E+00
 PARAMETER:  1.6903E-01 -1.0419E+00 -1.3632E+00  5.1733E-01 -1.2360E+00 -1.5098E-01  4.1834E-01 -1.2528E-01  6.0465E-01 -7.1190E-01
             9.0989E-01
 GRADIENT:   1.1676E+02  1.3749E+02  1.6762E+01  2.3934E+02 -8.2524E+01 -7.0326E+01 -4.4947E+00 -7.0306E+00  2.3064E-02 -7.8056E+00
             1.0770E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2974.88709525466        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      506
 NPARAMETR:  1.0247E+00  2.3460E-01  2.0506E-01  1.1768E+00  2.3732E-01  9.3823E-01  1.4161E+00  9.4775E-01  1.2788E+00  4.4230E-01
             2.0103E+00
 PARAMETER:  1.2442E-01 -1.3499E+00 -1.4845E+00  2.6278E-01 -1.3384E+00  3.6243E-02  4.4792E-01  4.6331E-02  3.4594E-01 -7.1577E-01
             7.9831E-01
 GRADIENT:  -2.2456E+01 -4.0095E+01  4.0017E+01  1.2081E+02  2.0225E+01  1.1047E+01  5.3117E+00 -6.3675E+00 -3.0058E+01 -4.2252E+00
            -5.2788E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2976.67507844059        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  1.0257E+00  2.4519E-01  2.0601E-01  1.1712E+00  2.3842E-01  8.8689E-01  1.3427E+00  9.5716E-01  1.2765E+00  4.4344E-01
             2.0183E+00
 PARAMETER:  1.2536E-01 -1.3057E+00 -1.4798E+00  2.5800E-01 -1.3337E+00 -2.0039E-02  3.9470E-01  5.6215E-02  3.4413E-01 -7.1319E-01
             8.0225E-01
 GRADIENT:   8.2227E+01  7.5577E+01  1.3329E+02  1.7091E+02  5.1943E+02 -3.8975E+00  1.3495E+00 -3.4926E+00 -7.2851E+00  6.6892E-01
            -3.3172E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2977.01320727669        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  1.0257E+00  2.4813E-01  2.0554E-01  1.1707E+00  2.3849E-01  8.9390E-01  1.3620E+00  9.5716E-01  1.2766E+00  4.4349E-01
             2.0215E+00
 PARAMETER:  1.2538E-01 -1.2938E+00 -1.4821E+00  2.5759E-01 -1.3334E+00 -1.2164E-02  4.0894E-01  5.6218E-02  3.4423E-01 -7.1309E-01
             8.0386E-01
 GRADIENT:  -2.2781E+01  2.9412E+00  3.3919E+01  1.1225E+02 -4.4616E+00 -7.0446E+00  4.8522E-01 -3.6562E+00 -2.8696E+01 -4.1098E+00
            -4.3992E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2984.21691555444        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1040            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0350E+00  2.4879E-01  2.0151E-01  1.0700E+00  2.3561E-01  9.1022E-01  1.3541E+00  9.7536E-01  1.3808E+00  4.8443E-01
             2.0615E+00
 PARAMETER:  1.3441E-01 -1.2912E+00 -1.5019E+00  1.6768E-01 -1.3456E+00  5.9327E-03  4.0313E-01  7.5048E-02  4.2269E-01 -6.2478E-01
             8.2345E-01
 GRADIENT:   1.1078E+02  9.7278E+01  1.3355E+02  3.9129E+01  5.3225E+02  6.7468E+00  7.3501E+00  3.1375E+00  2.7849E+01  7.2275E+00
             2.0503E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2984.35365253928        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1114
 NPARAMETR:  1.0380E+00  2.4682E-01  1.9893E-01  1.0356E+00  2.3236E-01  9.0845E-01  1.3525E+00  9.8301E-01  1.4318E+00  4.9566E-01
             2.0584E+00
 PARAMETER:  1.3733E-01 -1.2991E+00 -1.5148E+00  1.3493E-01 -1.3595E+00  3.9812E-03  4.0193E-01  8.2868E-02  4.5890E-01 -6.0187E-01
             8.2193E-01
 GRADIENT:   1.2154E+02  1.0282E+02  1.4626E+02  3.2936E+00  5.3695E+02  6.1268E+00  7.7671E+00  3.7237E+00  3.7547E+01  8.8991E+00
             2.1419E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2984.38718362534        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:     1223
 NPARAMETR:  1.0401E+00  2.4393E-01  1.9522E-01  1.0117E+00  2.2855E-01  9.0769E-01  1.3520E+00  9.9290E-01  1.4743E+00  5.0252E-01
             2.0545E+00
 PARAMETER:  1.3933E-01 -1.3109E+00 -1.5337E+00  1.1165E-01 -1.3760E+00  3.1440E-03  4.0156E-01  9.2879E-02  4.8817E-01 -5.8812E-01
             8.2005E-01
 GRADIENT:   1.6622E+01  2.4719E+01  5.9468E+01 -4.0114E+01 -7.0359E-01  8.9180E-02  3.8113E+00  3.5748E+00  2.0975E+01  4.5264E+00
             7.2952E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2986.47725377918        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  1.0464E+00  2.3633E-01  1.8173E-01  9.7694E-01  2.2316E-01  9.0725E-01  1.3278E+00  9.5513E-01  1.4927E+00  5.2129E-01
             2.0392E+00
 PARAMETER:  1.4539E-01 -1.3425E+00 -1.6052E+00  7.6665E-02 -1.3999E+00  2.6643E-03  3.8355E-01  5.4088E-02  5.0055E-01 -5.5146E-01
             8.1254E-01
 GRADIENT:   3.3020E+01  1.8290E+01  1.9718E+01 -5.8562E+01  5.8213E+01 -2.8372E-01  1.1725E+00  8.5033E-01  1.7534E+01  3.2382E+00
            -1.9834E+00

0ITERATION NO.:   51    OBJECTIVE VALUE:  -2986.47725377918        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1434
 NPARAMETR:  1.0468E+00  2.3567E-01  1.8233E-01  9.7673E-01  2.2379E-01  9.0764E-01  1.3228E+00  9.5493E-01  1.4911E+00  5.2188E-01
             2.0426E+00
 PARAMETER:  1.4539E-01 -1.3425E+00 -1.6052E+00  7.6665E-02 -1.3999E+00  2.6643E-03  3.8355E-01  5.4088E-02  5.0055E-01 -5.5146E-01
             8.1254E-01
 GRADIENT:  -4.2803E+04  2.3331E+03 -3.8477E+03  6.2216E+04 -4.3794E+03 -1.9711E-01  9.5581E-01  3.1134E+04  6.2366E+03 -1.1288E+04
            -7.7726E+03
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         1.9         1.6         2.3         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1434
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2384E-02  1.3040E-03 -2.5467E-03  3.7930E-02  1.3032E-03
 SE:             2.9554E-02  2.5032E-02  2.0872E-02  2.6813E-02  1.9210E-02
 N:                     100         100         100         100         100

 P VAL.:         6.7518E-01  9.5845E-01  9.0289E-01  1.5718E-01  9.4591E-01

 ETASHRINKSD(%)  9.9111E-01  1.6140E+01  3.0077E+01  1.0173E+01  3.5644E+01
 ETASHRINKVR(%)  1.9724E+00  2.9675E+01  5.1107E+01  1.9312E+01  5.8583E+01
 EBVSHRINKSD(%)  1.2574E+00  1.4576E+01  2.9058E+01  3.5776E+00  3.5730E+01
 EBVSHRINKVR(%)  2.4991E+00  2.7028E+01  4.9673E+01  7.0272E+00  5.8693E+01
 RELATIVEINF(%)  9.7482E+01  1.9219E+01  5.9402E+00  5.4626E+01  3.8145E+00
 EPSSHRINKSD(%)  2.1202E+01
 EPSSHRINKVR(%)  3.7909E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2986.4772537791837     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1332.3878940107729     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2986.477       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.36E-01  1.82E-01  9.77E-01  2.23E-01  9.07E-01  1.33E+00  9.55E-01  1.49E+00  5.21E-01  2.04E+00
 


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
 #CPUT: Total CPU Time in Seconds,      106.040
Stop Time:
Sat Oct 23 14:12:09 CDT 2021
