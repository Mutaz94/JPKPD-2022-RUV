Thu Sep 30 05:50:54 CDT 2021
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
$DATA ../../../../data/spa2/A2/dat62.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1281.91222419744        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9449E+02  1.2317E+02  1.1209E+02  7.7263E+01  2.1519E+02  6.6319E+01 -2.1806E+01 -1.6585E+02 -7.6489E+01 -1.0446E+02
            -2.0230E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1980.75034551008        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9699E-01  9.1356E-01  9.0232E-01  1.0177E+00  8.5222E-01  8.7116E-01  9.3096E-01  1.0373E+00  1.0003E+00  1.0013E+00
             2.1823E+00
 PARAMETER:  9.6989E-02  9.5957E-03 -2.7849E-03  1.1756E-01 -5.9914E-02 -3.7926E-02  2.8462E-02  1.3662E-01  1.0031E-01  1.0127E-01
             8.8036E-01
 GRADIENT:  -4.4733E+00  2.8664E+00  1.8496E+01 -2.5209E+01  1.6608E+01 -4.9941E+00  1.2597E+01  2.5707E+00  9.5286E-02  3.7853E+00
            -6.1809E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1990.24038909316        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      204
 NPARAMETR:  9.8756E-01  6.8590E-01  6.2665E-01  1.1579E+00  6.1233E-01  9.5349E-01  6.0246E-01  8.1315E-01  9.8572E-01  8.2761E-01
             2.1114E+00
 PARAMETER:  8.7484E-02 -2.7702E-01 -3.6737E-01  2.4658E-01 -3.9048E-01  5.2369E-02 -4.0674E-01 -1.0684E-01  8.5618E-02 -8.9214E-02
             8.4733E-01
 GRADIENT:  -1.0892E+02  9.7808E+00  1.0202E+01  4.1503E+01  2.1133E+01  1.5541E+01 -6.9114E-01  1.1321E+00  2.2468E+00 -1.8571E+01
            -9.2537E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2001.12122152722        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      384
 NPARAMETR:  1.0286E+00  6.5083E-01  5.4586E-01  1.1706E+00  5.7644E-01  9.2984E-01  3.2760E-01  2.9143E-01  9.6455E-01  9.5925E-01
             2.2611E+00
 PARAMETER:  1.2822E-01 -3.2951E-01 -5.0539E-01  2.5750E-01 -4.5088E-01  2.7262E-02 -1.0160E+00 -1.1330E+00  6.3905E-02  5.8399E-02
             9.1586E-01
 GRADIENT:  -1.9332E+01 -1.4449E+01 -2.6529E+01  6.1669E+01  7.2809E+01  1.2727E+01 -4.4760E-02  9.0593E-01 -4.6492E+00 -7.5412E+00
            -4.1319E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2007.75347656464        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      559
 NPARAMETR:  1.0342E+00  5.1745E-01  3.9321E-01  1.1428E+00  4.2176E-01  8.9762E-01  6.1341E-02  1.5868E-01  9.5684E-01  9.1423E-01
             2.2483E+00
 PARAMETER:  1.3363E-01 -5.5885E-01 -8.3340E-01  2.3345E-01 -7.6331E-01 -8.0118E-03 -2.6913E+00 -1.7409E+00  5.5883E-02  1.0332E-02
             9.1018E-01
 GRADIENT:  -6.3147E+00 -6.6104E-01 -2.1243E+00  1.5682E+00  3.5725E+00 -2.0330E+00 -9.1885E-03  3.6955E-01  7.9912E-01 -9.0656E-01
            -2.7613E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2007.95344129213        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  1.0367E+00  5.1590E-01  3.9178E-01  1.1415E+00  4.1962E-01  9.0246E-01  5.1399E-02  2.9897E-02  9.5258E-01  9.1781E-01
             2.2506E+00
 PARAMETER:  1.3605E-01 -5.6184E-01 -8.3706E-01  2.3233E-01 -7.6840E-01 -2.6343E-03 -2.8681E+00 -3.4100E+00  5.1419E-02  1.4236E-02
             9.1121E-01
 GRADIENT:   9.8636E-02  9.9632E-01  1.9624E-01 -6.7179E-01 -9.6145E-01  5.2699E-02 -5.8297E-03  1.3267E-02 -3.6911E-01 -3.0827E-01
             3.0164E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2007.96012438910        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  1.0368E+00  5.1528E-01  3.9201E-01  1.1420E+00  4.1953E-01  9.0235E-01  5.0127E-02  1.0000E-02  9.5353E-01  9.1879E-01
             2.2499E+00
 PARAMETER:  1.3618E-01 -5.6305E-01 -8.3647E-01  2.3279E-01 -7.6862E-01 -2.7511E-03 -2.8932E+00 -4.5315E+00  5.2420E-02  1.5299E-02
             9.1090E-01
 GRADIENT:   4.8027E-01  4.2300E-01  4.3357E-01 -2.9879E-01 -8.5443E-01  1.3532E-02 -5.2252E-03  4.0534E-04  8.5249E-03  7.1552E-03
            -5.1953E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2008.06318017470        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1099
 NPARAMETR:  1.0370E+00  5.1539E-01  3.9198E-01  1.1433E+00  4.1946E-01  9.0256E-01  4.2247E-01  1.0000E-02  9.5061E-01  8.9928E-01
             2.2439E+00
 PARAMETER:  1.3629E-01 -5.6284E-01 -8.3656E-01  2.3392E-01 -7.6878E-01 -2.5203E-03 -7.6163E-01 -5.0319E+00  4.9349E-02 -6.1617E-03
             9.0823E-01
 GRADIENT:   9.5570E-01 -4.6743E-01  5.9932E-01 -6.1848E-01  8.4053E-01 -2.2440E-02  1.0327E-01  0.0000E+00  1.6766E-01  2.1021E+00
             7.0995E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2008.07637403611        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1276
 NPARAMETR:  1.0370E+00  5.1521E-01  3.9071E-01  1.1431E+00  4.1836E-01  9.0280E-01  4.3622E-01  1.0000E-02  9.5006E-01  8.8842E-01
             2.2418E+00
 PARAMETER:  1.3631E-01 -5.6318E-01 -8.3978E-01  2.3375E-01 -7.7141E-01 -2.2489E-03 -7.2962E-01 -5.0629E+00  4.8771E-02 -1.8311E-02
             9.0729E-01
 GRADIENT:   1.0268E+00  1.2391E+00  1.4137E+00 -3.2611E-01 -9.2061E-01  2.1556E-02 -4.7031E-02  0.0000E+00 -7.5489E-02 -5.9496E-01
            -5.3047E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2008.08110523813        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1433
 NPARAMETR:  1.0360E+00  5.1400E-01  3.8955E-01  1.1428E+00  4.1762E-01  9.0263E-01  4.3820E-01  1.0000E-02  9.5008E-01  8.8925E-01
             2.2420E+00
 PARAMETER:  1.3534E-01 -5.6553E-01 -8.4275E-01  2.3352E-01 -7.7319E-01 -2.4434E-03 -7.2508E-01 -5.0629E+00  4.8790E-02 -1.7375E-02
             9.0736E-01
 GRADIENT:  -1.5426E+00  1.5744E-01  3.1523E-01 -6.1185E-01  9.4980E-01 -7.8421E-02 -1.4260E-02  0.0000E+00 -4.4143E-02 -6.6516E-02
            -1.7206E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -2008.08286743949        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     1524
 NPARAMETR:  1.0368E+00  5.1397E-01  3.8943E-01  1.1430E+00  4.1733E-01  9.0285E-01  4.4214E-01  1.0000E-02  9.5023E-01  8.8947E-01
             2.2422E+00
 PARAMETER:  1.3610E-01 -5.6560E-01 -8.4307E-01  2.3368E-01 -7.7388E-01 -2.2042E-03 -7.1614E-01 -5.0629E+00  4.8949E-02 -1.7131E-02
             9.0747E-01
 GRADIENT:   4.3583E-01  7.5871E-01  7.7156E-01 -2.6107E-01 -1.8776E-01  7.7113E-03  7.8696E-03  0.0000E+00  2.7516E-02  1.8901E-01
            -1.5017E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1524
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.3954E-04 -7.2382E-03 -4.7764E-05 -4.3042E-03 -3.4618E-03
 SE:             2.9404E-02  8.4216E-03  1.8344E-04  2.7737E-02  2.7195E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7993E-01  3.9007E-01  7.9457E-01  8.7668E-01  8.9871E-01

 ETASHRINKSD(%)  1.4923E+00  7.1787E+01  9.9385E+01  7.0776E+00  8.8919E+00
 ETASHRINKVR(%)  2.9624E+00  9.2040E+01  9.9996E+01  1.3654E+01  1.6993E+01
 EBVSHRINKSD(%)  1.7149E+00  7.2279E+01  9.9348E+01  6.7070E+00  8.9428E+00
 EBVSHRINKVR(%)  3.4004E+00  9.2316E+01  9.9996E+01  1.2964E+01  1.7086E+01
 RELATIVEINF(%)  9.6485E+01  1.3810E+00  6.6759E-04  7.2324E+01  6.2393E+00
 EPSSHRINKSD(%)  2.5228E+01
 EPSSHRINKVR(%)  4.4091E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2008.0828674394907     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -905.35662759388356     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2008.083       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  5.14E-01  3.89E-01  1.14E+00  4.17E-01  9.03E-01  4.42E-01  1.00E-02  9.50E-01  8.89E-01  2.24E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.04E+01  1.63E+03
 
 TH 3
+        1.41E+01  7.36E+02  2.76E+03
 
 TH 4
+       -3.38E+01  3.02E+02 -3.53E+02  8.21E+02
 
 TH 5
+        1.07E+01 -2.53E+03 -3.58E+03 -3.15E+01  6.70E+03
 
 TH 6
+        1.44E+00 -4.50E+00  1.88E+01 -7.43E+00 -7.08E-01  2.27E+02
 
 TH 7
+        2.58E-01 -7.54E+00 -2.94E+00 -4.70E+00  7.76E+00 -1.57E-01  5.44E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.97E+00 -5.28E+01  4.30E+01 -3.52E+00  2.53E+00 -2.63E-01  2.51E+00  0.00E+00  1.64E+02
 
 TH10
+       -8.97E-02 -1.82E+01 -6.05E+01  8.95E+00 -2.15E+01  1.09E-01  2.22E+01  0.00E+00  5.43E-01  1.75E+02
 
 TH11
+       -1.61E+01 -1.45E+01 -2.38E+01 -1.20E+01  2.65E+01  3.56E+00  4.98E+00  0.00E+00  9.64E+00  8.31E+00  1.21E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.132
Stop Time:
Thu Sep 30 05:51:28 CDT 2021
