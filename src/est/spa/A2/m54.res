Wed Sep 29 12:53:26 CDT 2021
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
$DATA ../../../../data/spa/A2/dat54.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -520.860282283630        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7108E+02  1.0523E+02  1.1565E+02  5.0399E+01  1.0064E+02  1.9967E+01 -1.7444E+01 -1.8257E+02 -1.4103E+02 -4.3634E+01
            -1.9090E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1288.79599383019        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9992E-01  9.2538E-01  9.4736E-01  1.0210E+00  9.4156E-01  8.8015E-01  9.5776E-01  1.1000E+00  1.1754E+00  7.8693E-01
             2.3615E+00
 PARAMETER:  9.9923E-02  2.2447E-02  4.5923E-02  1.2081E-01  3.9786E-02 -2.7665E-02  5.6839E-02  1.9532E-01  2.6162E-01 -1.3962E-01
             9.5930E-01
 GRADIENT:  -7.6317E+01 -2.1933E+00  6.6048E+00 -7.4786E+00  3.5233E+01 -7.3857E+01  2.5097E+00 -1.0245E+01 -8.0583E+00  5.1781E+00
            -2.4161E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1339.09722058069        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0101E+00  9.7236E-01  9.1120E-01  9.5157E-01  9.7429E-01  1.0325E+00  9.5434E-02  1.0758E+00  1.4533E+00  1.3816E-01
             3.1953E+00
 PARAMETER:  1.1006E-01  7.1974E-02  7.0106E-03  5.0363E-02  7.3955E-02  1.3197E-01 -2.2493E+00  1.7308E-01  4.7383E-01 -1.8794E+00
             1.2617E+00
 GRADIENT:  -5.2753E+01 -5.4920E+01 -9.9607E+00 -4.3669E+01  4.5694E+01  4.9935E+00  1.9691E-01 -3.0142E+00  2.4867E+01  3.8268E-01
            -7.9738E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1346.03082073518        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0403E+00  9.0085E-01  7.4447E-01  1.0399E+00  8.0956E-01  1.0076E+00  2.6971E-01  1.0539E+00  1.2286E+00  1.1649E-01
             3.2253E+00
 PARAMETER:  1.3951E-01 -4.4181E-03 -1.9508E-01  1.3916E-01 -1.1127E-01  1.0753E-01 -1.2104E+00  1.5250E-01  3.0587E-01 -2.0500E+00
             1.2710E+00
 GRADIENT:  -2.1888E-01  7.8725E+00 -3.8451E-01  7.0152E+00 -1.5773E+00 -2.4740E+00  1.0105E-01  2.9664E-02 -3.9529E-01  3.7516E-01
            -1.2388E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1347.56507753141        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0624E+00  7.8848E-01  7.8220E-01  1.1150E+00  7.8501E-01  1.0141E+00  2.1109E-01  9.8139E-01  1.1636E+00  9.1734E-02
             3.2842E+00
 PARAMETER:  1.6048E-01 -1.3765E-01 -1.4564E-01  2.0885E-01 -1.4205E-01  1.1403E-01 -1.4555E+00  8.1220E-02  2.5156E-01 -2.2889E+00
             1.2891E+00
 GRADIENT:  -1.9652E+01  2.0422E+00  1.2566E+00 -4.0632E+00 -3.0446E+00 -1.9361E+00  3.0356E-02  6.0682E-03 -5.1702E-01  2.2007E-01
            -8.3543E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1348.49712061512        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  1.0665E+00  4.6810E-01  7.5091E-01  1.3092E+00  6.5270E-01  1.0185E+00  2.6317E-02  8.3661E-01  9.7940E-01  1.9247E-02
             3.2716E+00
 PARAMETER:  1.6440E-01 -6.5907E-01 -1.8647E-01  3.6942E-01 -3.2664E-01  1.1835E-01 -3.5375E+00 -7.8399E-02  7.9182E-02 -3.8504E+00
             1.2853E+00
 GRADIENT:  -7.6139E+00  8.4801E+00  9.2036E+00  1.5872E+01 -1.5931E+01  7.5728E-01  4.0804E-04 -7.8882E-01 -2.8012E+00  1.0108E-02
            -3.5994E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1349.20013775293        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0658E+00  2.8203E-01  6.9609E-01  1.3874E+00  5.7931E-01  1.0085E+00  1.0000E-02  7.6846E-01  8.9996E-01  1.0000E-02
             3.2908E+00
 PARAMETER:  1.6369E-01 -1.1657E+00 -2.6227E-01  4.2740E-01 -4.4592E-01  1.0845E-01 -7.4869E+00 -1.6337E-01 -5.4015E-03 -6.1376E+00
             1.2911E+00
 GRADIENT:  -2.9892E-01  1.9289E+00  3.3745E+00  2.3292E+00 -4.8410E+00 -8.8526E-01  0.0000E+00 -4.9067E-01 -1.7537E+00  0.0000E+00
             1.6867E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1349.78777741803        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  1.0619E+00  1.0295E-01  6.1635E-01  1.4543E+00  4.9567E-01  1.0171E+00  1.0000E-02  8.4692E-01  8.4926E-01  1.0000E-02
             3.2023E+00
 PARAMETER:  1.6005E-01 -2.1735E+00 -3.8393E-01  4.7456E-01 -6.0184E-01  1.1691E-01 -1.7041E+01 -6.6145E-02 -6.3386E-02 -1.1165E+01
             1.2639E+00
 GRADIENT:   2.8000E+00  1.1338E+00  4.6793E+00  9.3573E+00 -9.3214E+00  8.8114E-01  0.0000E+00 -3.8200E-01 -2.0451E+00  0.0000E+00
            -3.8747E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1350.09064562639        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  1.0557E+00  2.3547E-02  6.2547E-01  1.4929E+00  4.9017E-01  1.0115E+00  1.0000E-02  8.8608E-01  8.2571E-01  1.0000E-02
             3.2058E+00
 PARAMETER:  1.5423E-01 -3.6487E+00 -3.6926E-01  5.0069E-01 -6.1299E-01  1.1144E-01 -3.2822E+01 -2.0944E-02 -9.1511E-02 -1.8962E+01
             1.2649E+00
 GRADIENT:  -1.4276E+00  1.4044E-01 -2.0775E-01  6.5901E+00 -7.1129E-01  1.0459E-01  0.0000E+00  1.9277E-02 -5.3115E-01  0.0000E+00
            -4.3014E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1350.15351683181        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  1.0555E+00  1.0000E-02  6.1218E-01  1.4900E+00  4.8073E-01  1.0110E+00  1.0000E-02  8.8141E-01  8.2467E-01  1.0000E-02
             3.2042E+00
 PARAMETER:  1.5398E-01 -4.6473E+00 -3.9072E-01  4.9879E-01 -6.3245E-01  1.1091E-01 -4.3744E+01 -2.6229E-02 -9.2773E-02 -2.4364E+01
             1.2645E+00
 GRADIENT:  -2.8053E-01  0.0000E+00 -3.4135E-01 -3.6789E-01  5.5383E-01 -7.8798E-03  0.0000E+00  2.2544E-02 -3.4540E-02  0.0000E+00
             9.1971E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1350.15351683181        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1253
 NPARAMETR:  1.0555E+00  1.0000E-02  6.1218E-01  1.4900E+00  4.8073E-01  1.0110E+00  1.0000E-02  8.8141E-01  8.2467E-01  1.0000E-02
             3.2042E+00
 PARAMETER:  1.5398E-01 -4.6473E+00 -3.9072E-01  4.9879E-01 -6.3245E-01  1.1091E-01 -4.3744E+01 -2.6229E-02 -9.2773E-02 -2.4364E+01
             1.2645E+00
 GRADIENT:  -2.8053E-01  0.0000E+00 -3.4135E-01 -3.6789E-01  5.5383E-01 -7.8798E-03  0.0000E+00  2.2544E-02 -3.4540E-02  0.0000E+00
             9.1971E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1253
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1968E-04 -2.6492E-06 -1.9031E-03 -9.0344E-03 -3.9921E-05
 SE:             2.8894E-02  1.4927E-06  1.6264E-02  2.6132E-02  2.3776E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9393E-01  7.5942E-02  9.0685E-01  7.2955E-01  8.6666E-01

 ETASHRINKSD(%)  3.2011E+00  9.9995E+01  4.5515E+01  1.2455E+01  9.9203E+01
 ETASHRINKVR(%)  6.2997E+00  1.0000E+02  7.0314E+01  2.3359E+01  9.9994E+01
 EBVSHRINKSD(%)  3.2459E+00  9.9995E+01  4.5428E+01  1.2227E+01  9.9125E+01
 EBVSHRINKVR(%)  6.3865E+00  1.0000E+02  7.0219E+01  2.2959E+01  9.9992E+01
 RELATIVEINF(%)  8.2725E+01  1.6588E-08  1.3737E+00  1.0690E+01  2.3801E-04
 EPSSHRINKSD(%)  2.7164E+01
 EPSSHRINKVR(%)  4.6950E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1350.1535168318146     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -615.00269026807644     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1350.154       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.00E-02  6.12E-01  1.49E+00  4.81E-01  1.01E+00  1.00E-02  8.81E-01  8.25E-01  1.00E-02  3.20E+00
 


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
+        9.15E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.74E+01  0.00E+00  1.13E+03
 
 TH 4
+       -6.36E+01  0.00E+00 -5.66E+01  5.66E+02
 
 TH 5
+        9.50E+01  0.00E+00 -1.88E+03 -2.82E+02  3.53E+03
 
 TH 6
+       -1.90E+00  0.00E+00  7.18E+00 -1.31E+01  1.03E+01  1.69E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -5.22E+00  0.00E+00 -4.91E+01 -2.42E+00  5.91E+01 -2.03E+00  0.00E+00  2.19E+01
 
 TH 9
+        6.82E+00  0.00E+00  6.95E+00 -1.88E+01  4.54E+01 -2.12E+00  0.00E+00  9.83E-01  1.71E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.26E+01  0.00E+00 -2.60E+00 -6.68E+00 -1.09E+01  4.71E+00  0.00E+00  1.41E+01  1.26E+01  0.00E+00  3.20E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.702
Stop Time:
Wed Sep 29 12:53:49 CDT 2021
