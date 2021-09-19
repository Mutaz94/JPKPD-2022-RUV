Sat Sep 18 09:06:25 CDT 2021
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
$DATA ../../../../data/spa/A1/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1236.24187954394        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.1137E+01 -2.9777E+01  5.7126E+01 -1.0483E+02  1.9309E+01  3.7570E+00 -3.6238E+01 -1.5764E+01 -3.9253E+01 -7.6563E+01
            -7.3173E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1461.22219245475        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0388E+00  1.0265E+00  9.6334E-01  1.1002E+00  1.0392E+00  9.5173E-01  1.1567E+00  9.4446E-01  1.0267E+00  1.2515E+00
             2.0280E+00
 PARAMETER:  1.3809E-01  1.2616E-01  6.2651E-02  1.9547E-01  1.3845E-01  5.0530E-02  2.4561E-01  4.2856E-02  1.2634E-01  3.2433E-01
             8.0705E-01
 GRADIENT:   4.6227E-01  1.0462E+01 -2.2394E+01  4.4548E+01  2.3090E+01 -8.7488E+00 -5.6866E+00  5.0220E+00  1.9675E+00 -5.6356E+00
            -1.0463E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1463.11491313932        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0376E+00  8.9607E-01  1.1629E+00  1.1981E+00  1.0643E+00  9.8458E-01  1.3699E+00  6.6010E-01  9.6520E-01  1.3921E+00
             2.0301E+00
 PARAMETER:  1.3696E-01 -9.7330E-03  2.5091E-01  2.8076E-01  1.6228E-01  8.4457E-02  4.1473E-01 -3.1536E-01  6.4579E-02  4.3083E-01
             8.0811E-01
 GRADIENT:   5.6306E+00  2.6388E+01  5.9726E+00  5.0954E+01 -2.6708E+00  6.8060E+00 -5.7480E-01 -4.9772E-01  1.9648E+00 -1.9947E-01
            -1.6983E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1466.65431267361        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0330E+00  5.6838E-01  6.6510E-01  1.2963E+00  6.2228E-01  9.6302E-01  1.9997E+00  8.4758E-02  7.9287E-01  9.0277E-01
             2.0207E+00
 PARAMETER:  1.3251E-01 -4.6497E-01 -3.0782E-01  3.5951E-01 -3.7436E-01  6.2315E-02  7.9301E-01 -2.3680E+00 -1.3209E-01 -2.2849E-03
             8.0346E-01
 GRADIENT:  -1.2124E+01  1.1823E+01  7.9544E+00  2.4100E+01 -1.0793E+01 -3.7978E+00 -6.9116E-01  7.6552E-02 -9.5911E+00 -1.2406E+00
            -6.8737E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1467.29627360743        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0412E+00  4.8740E-01  5.0748E-01  1.2941E+00  4.9608E-01  9.7886E-01  2.1394E+00  2.7163E-02  8.1967E-01  7.4449E-01
             2.0176E+00
 PARAMETER:  1.4036E-01 -6.1866E-01 -5.7829E-01  3.5782E-01 -6.0103E-01  7.8628E-02  8.6053E-01 -3.5059E+00 -9.8850E-02 -1.9506E-01
             8.0191E-01
 GRADIENT:   2.2937E+00  9.2558E+00  2.5530E+00  3.3009E+01 -4.3641E+00  1.1495E+00  7.0022E-01  1.0544E-02  8.6938E-01 -1.1896E+00
            -2.4997E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1467.81543906915        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      485
 NPARAMETR:  1.0350E+00  4.1667E-01  6.7004E-01  1.3647E+00  5.8829E-01  9.6113E-01  2.4630E+00  3.6259E-02  7.9246E-01  9.1645E-01
             2.0246E+00
 PARAMETER:  1.3440E-01 -7.7546E-01 -3.0042E-01  4.1097E-01 -4.3053E-01  6.0352E-02  1.0014E+00 -3.2171E+00 -1.3261E-01  1.2749E-02
             8.0538E-01
 GRADIENT:  -1.3304E+01  3.0696E+00  5.0243E+00  9.0151E-01 -7.6471E+00 -3.8400E+00  1.1832E-01  1.4032E-02 -4.8870E+00  2.4280E-01
            -4.2184E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1468.06499286168        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.0387E+00  3.4385E-01  6.9541E-01  1.4057E+00  5.9138E-01  9.6816E-01  2.7663E+00  2.2038E-02  8.0015E-01  9.4196E-01
             2.0404E+00
 PARAMETER:  1.3799E-01 -9.6754E-01 -2.6326E-01  4.4055E-01 -4.2530E-01  6.7647E-02  1.1175E+00 -3.7150E+00 -1.2295E-01  4.0211E-02
             8.1315E-01
 GRADIENT:   4.7843E-02 -1.3108E-01 -3.6634E-01 -2.1159E-01  5.1466E-01 -6.8009E-03 -2.9924E-02  4.8327E-03 -4.8176E-03 -4.1671E-02
             2.3456E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1468.06519800498        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      837
 NPARAMETR:  1.0387E+00  3.4461E-01  6.9628E-01  1.4057E+00  5.9192E-01  9.6817E-01  2.7643E+00  2.1906E-02  8.0018E-01  9.4296E-01
             2.0403E+00
 PARAMETER:  1.3797E-01 -9.6535E-01 -2.6201E-01  4.4054E-01 -4.2439E-01  6.7656E-02  1.1168E+00 -3.7210E+00 -1.2292E-01  4.1265E-02
             8.1309E-01
 GRADIENT:   4.9805E-03 -1.0819E-02 -5.7907E-02 -4.3690E-02  6.3507E-02 -5.4924E-03  1.8572E-02  4.7616E-03 -2.3181E-02 -3.2267E-03
            -2.5557E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1468.06696606708        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1016
 NPARAMETR:  1.0388E+00  3.4771E-01  6.9705E-01  1.4045E+00  5.9314E-01  9.6823E-01  2.7495E+00  1.0000E-02  8.0062E-01  9.4353E-01
             2.0405E+00
 PARAMETER:  1.3805E-01 -9.5639E-01 -2.6090E-01  4.3971E-01 -4.2232E-01  6.7716E-02  1.1114E+00 -4.6713E+00 -1.2237E-01  4.1876E-02
             8.1322E-01
 GRADIENT:   3.1089E-03  3.8245E-02 -1.5336E-01  1.9193E-01  1.7899E-01 -1.0435E-02  2.9459E-03  0.0000E+00 -3.4946E-03  1.3517E-03
             1.8365E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1468.06705486558        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1192
 NPARAMETR:  1.0388E+00  3.4678E-01  6.9697E-01  1.4049E+00  5.9284E-01  9.6822E-01  2.7536E+00  1.0032E-02  8.0051E-01  9.4343E-01
             2.0404E+00
 PARAMETER:  1.3802E-01 -9.5907E-01 -2.6102E-01  4.3996E-01 -4.2284E-01  6.7705E-02  1.1129E+00 -4.5019E+00 -1.2251E-01  4.1772E-02
             8.1316E-01
 GRADIENT:   4.1290E-03  1.5768E-02 -4.1009E-02  5.4049E-02  6.2036E-02 -4.6733E-03 -8.4218E-03  8.7236E-04 -4.4670E-03 -1.0395E-02
            -1.0898E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1468.06706210442        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1286
 NPARAMETR:  1.0388E+00  3.4672E-01  6.9700E-01  1.4049E+00  5.9281E-01  9.6823E-01  2.7538E+00  1.0000E-02  8.0052E-01  9.4350E-01
             2.0405E+00
 PARAMETER:  1.3802E-01 -9.5923E-01 -2.6098E-01  4.3993E-01 -4.2288E-01  6.7715E-02  1.1130E+00 -4.5157E+00 -1.2249E-01  4.1841E-02
             8.1319E-01
 GRADIENT:   1.2948E-02  1.0261E-02  3.7541E-02 -8.8633E-02 -2.7017E-02  1.6636E-03 -6.7040E-03  0.0000E+00  1.5962E-03 -2.9744E-04
             4.2094E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1286
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9740E-04  3.4676E-02 -2.7291E-04 -3.0442E-02 -2.1580E-03
 SE:             2.9369E-02  1.7820E-02  1.6856E-04  2.4134E-02  2.0918E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7291E-01  5.1659E-02  1.0543E-01  2.0717E-01  9.1783E-01

 ETASHRINKSD(%)  1.6113E+00  4.0302E+01  9.9435E+01  1.9148E+01  2.9923E+01
 ETASHRINKVR(%)  3.1966E+00  6.4361E+01  9.9997E+01  3.4630E+01  5.0892E+01
 EBVSHRINKSD(%)  1.7533E+00  4.7119E+01  9.9381E+01  1.6572E+01  2.5341E+01
 EBVSHRINKVR(%)  3.4759E+00  7.2036E+01  9.9996E+01  3.0399E+01  4.4261E+01
 RELATIVEINF(%)  9.5185E+01  6.3722E+00  2.4590E-04  1.9921E+01  3.3199E+00
 EPSSHRINKSD(%)  3.6886E+01
 EPSSHRINKVR(%)  6.0166E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1468.0670621044212     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -732.91623554068303     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1468.067       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.47E-01  6.97E-01  1.40E+00  5.93E-01  9.68E-01  2.75E+00  1.00E-02  8.01E-01  9.43E-01  2.04E+00
 


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
+        1.06E+03
 
 TH 2
+       -3.73E+01  4.00E+02
 
 TH 3
+        2.47E+01  1.74E+02  9.43E+02
 
 TH 4
+       -2.88E+01  2.57E+02 -1.87E+02  6.35E+02
 
 TH 5
+        5.87E+00 -4.00E+02 -1.30E+03  8.11E+01  2.01E+03
 
 TH 6
+       -2.68E+00 -9.06E+00  7.93E+00 -7.97E+00 -2.65E+00  2.01E+02
 
 TH 7
+        1.62E+00  3.43E+01 -3.39E+00 -6.40E+00  5.92E-01  2.72E-01  7.91E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.81E+00 -2.30E+01  1.91E-01 -8.03E+00  1.35E+01 -4.58E+00  1.65E+00  0.00E+00  1.71E+02
 
 TH10
+       -8.27E+00  9.32E+00 -4.67E+01 -2.06E+01 -2.60E+01 -6.64E+00  1.69E+00  0.00E+00 -3.78E+00  7.49E+01
 
 TH11
+       -1.17E+01 -5.00E+00 -2.34E+01 -1.24E+01  1.20E+01  4.45E+00  8.18E-01  0.00E+00  1.26E+01  1.45E+01  6.14E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.383
Stop Time:
Sat Sep 18 09:06:49 CDT 2021
