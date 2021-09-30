Wed Sep 29 18:09:28 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1695.60969266850        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5183E+02 -6.4479E+01 -2.6974E+01 -2.7532E+01  8.0712E+01  4.8134E+01  7.6328E+00  9.5916E-01  3.0504E+01 -7.2440E+00
             5.4667E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1704.98409302048        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0368E+00  1.0679E+00  9.7833E-01  1.0432E+00  9.5465E-01  9.4742E-01  9.2600E-01  1.0129E+00  8.0838E-01  1.0154E+00
             1.0129E+00
 PARAMETER:  1.3618E-01  1.6570E-01  7.8088E-02  1.4231E-01  5.3594E-02  4.5986E-02  2.3115E-02  1.1281E-01 -1.1272E-01  1.1524E-01
             1.1278E-01
 GRADIENT:   1.1183E+01  1.5132E+01 -3.1332E-01  1.7603E+01  2.4352E+00 -1.0837E+01 -3.8457E+00 -3.0622E-01 -1.2466E+01  1.6958E+00
             1.0737E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1705.37094501556        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0395E+00  1.0190E+00  8.9787E-01  1.0655E+00  8.9248E-01  9.5894E-01  1.0037E+00  9.2622E-01  7.8433E-01  9.3579E-01
             1.0057E+00
 PARAMETER:  1.3874E-01  1.1879E-01 -7.7300E-03  1.6348E-01 -1.3748E-02  5.8072E-02  1.0369E-01  2.3355E-02 -1.4293E-01  3.3631E-02
             1.0570E-01
 GRADIENT:   1.6091E+01  9.8325E+00 -3.6662E+00  1.5768E+01  5.6919E+00 -6.1551E+00 -2.3040E+00  1.0095E+00 -1.0169E+01 -5.7978E-01
            -9.0843E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1705.90057990623        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0317E+00  9.3780E-01  8.9095E-01  1.1067E+00  8.5185E-01  9.7428E-01  1.0626E+00  8.2451E-01  8.0390E-01  9.2175E-01
             1.0037E+00
 PARAMETER:  1.3117E-01  3.5785E-02 -1.5468E-02  2.0139E-01 -6.0339E-02  7.3947E-02  1.6075E-01 -9.2961E-02 -1.1828E-01  1.8520E-02
             1.0374E-01
 GRADIENT:  -1.3032E+00  2.7920E+00  1.5956E+00  1.8655E+00 -2.2899E+00  4.0059E-01  2.0385E-01 -2.7976E-01  1.0814E-01 -2.6639E-01
             1.7156E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1705.90260400217        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0313E+00  8.8913E-01  9.2561E-01  1.1394E+00  8.4864E-01  9.7326E-01  1.1016E+00  8.4511E-01  7.8857E-01  9.2680E-01
             1.0040E+00
 PARAMETER:  1.3085E-01 -1.7509E-02  2.2693E-02  2.3050E-01 -6.4117E-02  7.2895E-02  1.9676E-01 -6.8294E-02 -1.3753E-01  2.3984E-02
             1.0398E-01
 GRADIENT:  -7.5769E-01  4.1324E+00  2.6857E+00  3.2089E+00 -4.2395E+00  2.3154E-01  1.1299E-01 -4.3575E-01 -5.1100E-02 -3.1616E-01
             1.3801E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1705.90372670117        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0312E+00  8.5001E-01  9.5327E-01  1.1655E+00  8.4629E-01  9.7206E-01  1.1335E+00  8.6230E-01  7.7690E-01  9.3141E-01
             1.0042E+00
 PARAMETER:  1.3075E-01 -6.2504E-02  5.2142E-02  2.5318E-01 -6.6899E-02  7.1659E-02  2.2534E-01 -4.8155E-02 -1.5244E-01  2.8941E-02
             1.0421E-01
 GRADIENT:   1.5550E-01  4.9358E+00  3.4296E+00  4.2970E+00 -5.7327E+00 -4.5488E-02 -2.9940E-02 -5.2364E-01 -2.5146E-01 -2.9833E-01
             2.4943E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1705.90836663216        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1075
 NPARAMETR:  1.0311E+00  7.9677E-01  9.9458E-01  1.2010E+00  8.4550E-01  9.7014E-01  1.1771E+00  8.9138E-01  7.6246E-01  9.4018E-01
             1.0045E+00
 PARAMETER:  1.3063E-01 -1.2719E-01  9.4564E-02  2.8315E-01 -6.7829E-02  6.9689E-02  2.6304E-01 -1.4989E-02 -1.7121E-01  3.8317E-02
             1.0450E-01
 GRADIENT:   1.6303E+00  5.4437E+00  4.0838E+00  5.3674E+00 -7.2822E+00 -5.0841E-01 -2.5381E-01 -5.7017E-01 -5.5078E-01 -1.9729E-01
             3.8926E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1705.92478895566        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  1.0309E+00  7.0305E-01  1.0767E+00  1.2634E+00  8.4951E-01  9.6625E-01  1.2518E+00  9.5425E-01  7.4114E-01  9.6123E-01
             1.0052E+00
 PARAMETER:  1.3045E-01 -2.5233E-01  1.7394E-01  3.3384E-01 -6.3100E-02  6.5665E-02  3.2457E-01  5.3169E-02 -1.9957E-01  6.0462E-02
             1.0517E-01
 GRADIENT:   4.7652E+00  5.5176E+00  4.6571E+00  6.5841E+00 -9.2083E+00 -1.4885E+00 -6.5294E-01 -5.4816E-01 -1.0649E+00  9.8184E-02
             6.4746E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1705.95894479224        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1425
 NPARAMETR:  1.0297E+00  5.8499E-01  1.2125E+00  1.3445E+00  8.6801E-01  9.6195E-01  1.3341E+00  1.0693E+00  7.2154E-01  9.9970E-01
             1.0061E+00
 PARAMETER:  1.2927E-01 -4.3616E-01  2.9266E-01  3.9605E-01 -4.1549E-02  6.1209E-02  3.8822E-01  1.6703E-01 -2.2637E-01  9.9699E-02
             1.0611E-01
 GRADIENT:   7.2746E+00  6.0195E+00  4.9646E+00  1.0128E+01 -1.1071E+01 -2.3954E+00 -7.7452E-01 -3.7841E-01 -1.3451E+00  6.6596E-01
             7.9908E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1706.00552485000        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1600
 NPARAMETR:  1.0272E+00  4.5780E-01  1.3856E+00  1.4331E+00  8.9619E-01  9.5914E-01  1.4093E+00  1.2275E+00  7.0430E-01  1.0441E+00
             1.0067E+00
 PARAMETER:  1.2680E-01 -6.8132E-01  4.2612E-01  4.5984E-01 -9.6017E-03  5.8286E-02  4.4309E-01  3.0494E-01 -2.5055E-01  1.4311E-01
             1.0665E-01
 GRADIENT:   7.1587E+00  5.9554E+00  3.5978E+00  1.5627E+01 -1.0208E+01 -2.6190E+00 -6.2262E-01  1.9399E-01 -1.0093E+00  1.4038E+00
             6.2152E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1706.06253850125        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1776
 NPARAMETR:  1.0239E+00  3.3916E-01  1.5377E+00  1.5128E+00  9.1679E-01  9.5917E-01  1.5047E+00  1.3650E+00  6.8581E-01  1.0703E+00
             1.0063E+00
 PARAMETER:  1.2361E-01 -9.8129E-01  5.3030E-01  5.1394E-01  1.3125E-02  5.8317E-02  5.0859E-01  4.1118E-01 -2.7715E-01  1.6796E-01
             1.0625E-01
 GRADIENT:   4.4166E+00  4.3816E+00  1.6778E+00  1.6264E+01 -5.9693E+00 -1.7714E+00 -5.1894E-01  1.4891E-02 -7.3109E-01  1.2033E+00
             7.2607E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1706.11911242678        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1958
 NPARAMETR:  1.0223E+00  2.7364E-01  1.6128E+00  1.5401E+00  9.3079E-01  9.6254E-01  2.0651E+00  1.4386E+00  6.7795E-01  1.0719E+00
             1.0065E+00
 PARAMETER:  1.2205E-01 -1.1959E+00  5.7798E-01  5.3183E-01  2.8283E-02  6.1823E-02  8.2518E-01  4.6364E-01 -2.8868E-01  1.6939E-01
             1.0650E-01
 GRADIENT:   3.7291E+00 -1.7654E+00 -1.3808E+00 -2.8829E+01  6.2644E+00  1.7631E-01  8.2528E-01  1.2911E-01  7.0849E+00  2.0672E-01
             7.0599E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1706.25729809474        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2133
 NPARAMETR:  1.0203E+00  2.5840E-01  1.6158E+00  1.5595E+00  9.2247E-01  9.6200E-01  2.1058E+00  1.4420E+00  6.5764E-01  1.0639E+00
             1.0061E+00
 PARAMETER:  1.2011E-01 -1.2532E+00  5.7983E-01  5.4437E-01  1.9301E-02  6.1259E-02  8.4471E-01  4.6600E-01 -3.1910E-01  1.6196E-01
             1.0610E-01
 GRADIENT:  -6.8835E-01  1.6411E+00 -1.0637E-01 -1.2063E+00  2.9513E-01 -6.2855E-03 -2.7280E-02 -1.6035E-01 -8.5709E-01 -7.9877E-03
            -1.3409E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1706.28058617931        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2317
 NPARAMETR:  1.0204E+00  2.5037E-01  1.6151E+00  1.5599E+00  9.2099E-01  9.6177E-01  2.1124E+00  1.4439E+00  6.5899E-01  1.0631E+00
             1.0060E+00
 PARAMETER:  1.2017E-01 -1.2848E+00  5.7938E-01  5.4461E-01  1.7691E-02  6.1021E-02  8.4783E-01  4.6738E-01 -3.1705E-01  1.6116E-01
             1.0600E-01
 GRADIENT:  -1.3585E-01  2.5215E-01 -3.3915E-01 -1.3673E+01  1.5736E+00 -2.6302E-02 -5.0379E-02 -1.1057E-02 -7.3326E-02 -3.6440E-02
             4.6421E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1706.28362213421        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2497
 NPARAMETR:  1.0213E+00  2.4766E-01  1.6132E+00  1.5605E+00  9.1960E-01  9.6201E-01  2.1366E+00  1.4443E+00  6.5871E-01  1.0604E+00
             1.0059E+00
 PARAMETER:  1.2112E-01 -1.2957E+00  5.7821E-01  5.4499E-01  1.6185E-02  6.1269E-02  8.5923E-01  4.6765E-01 -3.1747E-01  1.5864E-01
             1.0586E-01
 GRADIENT:   2.2880E+00 -5.3248E-03 -3.4237E-01 -1.6342E+01  1.7179E+00  9.2723E-02 -9.4747E-03  6.5387E-02  1.0795E-01 -1.8269E-01
             1.3979E-02

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1706.28548966744        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     2594
 NPARAMETR:  1.0213E+00  2.4833E-01  1.6138E+00  1.5611E+00  9.1838E-01  9.6199E-01  2.1464E+00  1.4440E+00  6.5872E-01  1.0615E+00
             1.0059E+00
 PARAMETER:  1.2110E-01 -1.2930E+00  5.7861E-01  5.4537E-01  1.4856E-02  6.1245E-02  8.6380E-01  4.6741E-01 -3.1746E-01  1.5966E-01
             1.0585E-01
 GRADIENT:   1.0283E-02  1.7853E-01  5.2942E-01 -4.1698E-01 -3.4964E-01  2.7909E-03  1.3856E-03  4.7263E-02  1.1313E-01  1.0061E-01
             1.5080E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2594
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8254E-04  1.8323E-03 -3.6786E-02 -8.3003E-03 -4.1482E-02
 SE:             2.9841E-02  9.6346E-03  1.8197E-02  2.7562E-02  2.0675E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9245E-01  8.4917E-01  4.3226E-02  7.6330E-01  4.4818E-02

 ETASHRINKSD(%)  2.9819E-02  6.7723E+01  3.9037E+01  7.6629E+00  3.0736E+01
 ETASHRINKVR(%)  5.9629E-02  8.9582E+01  6.2836E+01  1.4739E+01  5.2025E+01
 EBVSHRINKSD(%)  4.4886E-01  6.8898E+01  4.2767E+01  7.7942E+00  2.6917E+01
 EBVSHRINKVR(%)  8.9570E-01  9.0326E+01  6.7244E+01  1.4981E+01  4.6588E+01
 RELATIVEINF(%)  9.5508E+01  2.5062E-01  8.0368E+00  2.5160E+00  7.5199E+00
 EPSSHRINKSD(%)  4.5225E+01
 EPSSHRINKVR(%)  6.9997E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1706.2854896674398     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -971.13466310370166     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.30
 Elapsed covariance  time in seconds:     6.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1706.285       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.48E-01  1.61E+00  1.56E+00  9.18E-01  9.62E-01  2.15E+00  1.44E+00  6.59E-01  1.06E+00  1.01E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.84E-02  1.40E+00  2.36E+00  9.53E-01  3.41E-01  7.06E-02  5.95E+00  2.31E+00  2.45E-01  2.86E-01  8.16E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.47E-03
 
 TH 2
+        3.41E-02  1.96E+00
 
 TH 3
+       -5.65E-02 -3.24E+00  5.58E+00
 
 TH 4
+       -2.32E-02 -1.33E+00  2.22E+00  9.07E-01
 
 TH 5
+       -7.47E-03 -4.32E-01  7.79E-01  2.97E-01  1.16E-01
 
 TH 6
+       -1.79E-04  2.29E-02 -4.11E-02 -1.58E-02 -5.97E-03  4.99E-03
 
 TH 7
+       -1.44E-01 -8.32E+00  1.38E+01  5.66E+00  1.85E+00 -9.86E-02  3.54E+01
 
 TH 8
+       -5.53E-02 -3.19E+00  5.43E+00  2.17E+00  7.48E-01 -4.01E-02  1.36E+01  5.33E+00
 
 TH 9
+        5.86E-03  3.37E-01 -5.63E-01 -2.30E-01 -7.55E-02  3.33E-03 -1.43E+00 -5.52E-01  6.02E-02
 
 TH10
+       -5.22E-03 -3.23E-01  5.80E-01  2.22E-01  8.73E-02 -4.52E-03  1.39E+00  5.51E-01 -5.73E-02  8.20E-02
 
 TH11
+       -1.23E-04 -1.72E-02  3.42E-02  1.25E-02  5.30E-03 -1.14E-03  7.48E-02  3.11E-02 -3.29E-03  3.21E-03  6.65E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.84E-02
 
 TH 2
+        6.35E-01  1.40E+00
 
 TH 3
+       -6.23E-01 -9.80E-01  2.36E+00
 
 TH 4
+       -6.35E-01 -9.99E-01  9.84E-01  9.53E-01
 
 TH 5
+       -5.72E-01 -9.06E-01  9.68E-01  9.15E-01  3.41E-01
 
 TH 6
+       -6.59E-02  2.31E-01 -2.47E-01 -2.35E-01 -2.48E-01  7.06E-02
 
 TH 7
+       -6.29E-01 -9.98E-01  9.83E-01  9.98E-01  9.14E-01 -2.35E-01  5.95E+00
 
 TH 8
+       -6.24E-01 -9.86E-01  9.96E-01  9.89E-01  9.52E-01 -2.46E-01  9.90E-01  2.31E+00
 
 TH 9
+        6.22E-01  9.82E-01 -9.70E-01 -9.82E-01 -9.03E-01  1.92E-01 -9.81E-01 -9.74E-01  2.45E-01
 
 TH10
+       -4.75E-01 -8.06E-01  8.58E-01  8.14E-01  8.96E-01 -2.24E-01  8.15E-01  8.34E-01 -8.16E-01  2.86E-01
 
 TH11
+       -3.93E-02 -1.51E-01  1.78E-01  1.60E-01  1.91E-01 -1.97E-01  1.54E-01  1.65E-01 -1.64E-01  1.37E-01  8.16E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.28E+03
 
 TH 2
+       -6.10E+01  5.92E+02
 
 TH 3
+        1.58E+01  4.54E+01  7.92E+01
 
 TH 4
+        1.25E+01  6.09E+02 -1.25E+01  9.38E+02
 
 TH 5
+       -3.18E+00 -1.89E+02 -1.81E+02 -4.33E+01  6.97E+02
 
 TH 6
+        1.61E+02 -9.12E+00  5.20E+00  7.52E+00 -2.19E+01  2.55E+02
 
 TH 7
+       -1.20E+01  4.50E+01  4.06E+00  1.23E+01 -2.37E+00 -1.88E+00  1.06E+01
 
 TH 8
+       -3.75E+00 -3.49E+01 -3.05E+01 -2.92E+01 -2.61E+00  6.26E+00 -8.85E+00  4.36E+01
 
 TH 9
+        5.58E+01 -6.15E+01  2.89E+01  2.43E+00 -9.45E+01  9.37E+01 -1.51E+00  3.87E+00  5.49E+02
 
 TH10
+       -2.41E+01 -9.09E+00  2.99E+00 -1.40E+01 -1.02E+02  8.47E+00 -4.59E+00  1.90E+01  3.84E+01  7.72E+01
 
 TH11
+       -1.73E+01 -6.22E+01 -1.17E+01 -7.48E+01  2.70E+00  3.36E+01 -2.78E+00  1.28E+01  3.15E+01  1.43E+01  1.73E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       39.644
Stop Time:
Wed Sep 29 18:10:09 CDT 2021
