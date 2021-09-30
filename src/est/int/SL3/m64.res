Wed Sep 29 04:34:55 CDT 2021
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
$DATA ../../../../data/int/SL3/dat64.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      984
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

 TOT. NO. OF OBS RECS:      884
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   271.230689817861        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6096E+02  5.5854E+00 -1.1005E+01  2.9534E+02  2.4193E+02  3.2094E+01 -1.3863E+02 -2.0889E+02 -2.3176E+02 -3.4886E+01
            -7.5796E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2334.73254290352        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0021E+00  1.4426E+00  1.0049E+00  8.3555E-01  1.1205E+00  9.9312E-01  1.0749E+00  1.0623E+00  1.0447E+00  1.0455E+00
             5.2118E+00
 PARAMETER:  1.0212E-01  4.6648E-01  1.0484E-01 -7.9665E-02  2.1376E-01  9.3092E-02  1.7222E-01  1.6046E-01  1.4372E-01  1.4453E-01
             1.7509E+00
 GRADIENT:  -1.3482E+02  2.0439E+01 -9.8561E+00  5.0363E+00  7.5368E+00 -2.5836E+01  3.0308E+01  3.4025E+00  1.6369E+01  2.7721E+00
             8.1722E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2420.84566241970        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.8944E-01  2.3756E+00  1.8829E+00  5.5855E-01  1.9234E+00  9.8412E-01  1.0022E+00  3.3702E+00  3.4442E+00  1.7040E+00
             4.2604E+00
 PARAMETER:  8.9381E-02  9.6525E-01  7.3279E-01 -4.8241E-01  7.5407E-01  8.3993E-02  1.0225E-01  1.3150E+00  1.3367E+00  6.3297E-01
             1.5494E+00
 GRADIENT:  -1.2624E+02  2.9390E+02 -7.9992E+00  9.2278E+01  5.8397E+01 -3.7335E+01  4.5425E+01  6.2317E+00  4.1096E+01 -5.1215E+00
             6.4427E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2624.31787003056        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.7159E-01  1.3939E+00  7.2800E+00  7.9779E-01  1.8170E+00  1.0375E+00  7.9419E-01  2.6701E-01  9.3212E-01  2.0780E+00
             2.7731E+00
 PARAMETER:  7.1176E-02  4.3210E-01  2.0851E+00 -1.2591E-01  6.9717E-01  1.3678E-01 -1.3044E-01 -1.2205E+00  2.9705E-02  8.3141E-01
             1.1199E+00
 GRADIENT:  -4.9988E+01 -2.1996E+01 -1.2141E+01 -3.6516E+01  4.9723E+01 -2.2008E+00 -1.0446E+01 -1.0298E-01 -1.7096E+01  2.6332E+01
            -5.6387E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2639.36561772255        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  9.9612E-01  1.0404E+00  5.2428E+01  1.0926E+00  1.8387E+00  1.0364E+00  8.0104E-01  1.9437E-02  9.0923E-01  2.0032E+00
             2.8002E+00
 PARAMETER:  9.6109E-02  1.3957E-01  4.0594E+00  1.8854E-01  7.0908E-01  1.3575E-01 -1.2185E-01 -3.8406E+00  4.8404E-03  7.9474E-01
             1.1297E+00
 GRADIENT:  -4.6570E+00  4.2679E+00 -2.7357E+00  2.3540E+01  2.3445E+01 -4.9925E-01 -5.2157E-01 -3.4258E-05 -5.5721E+00  6.2916E+00
            -2.0824E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2639.51942900491        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  9.9787E-01  1.0372E+00  5.8381E+01  1.0815E+00  1.7929E+00  1.0353E+00  7.7469E-01  1.6310E-02  9.4365E-01  1.9667E+00
             2.8163E+00
 PARAMETER:  9.7864E-02  1.3655E-01  4.1670E+00  1.7831E-01  6.8385E-01  1.3474E-01 -1.5530E-01 -4.0160E+00  4.2000E-02  7.7636E-01
             1.1354E+00
 GRADIENT:  -1.3068E+00 -8.7520E+00 -2.0945E+00  3.3382E+00  1.0268E+01 -6.0348E-01  1.3557E+00 -1.9240E-05 -7.6185E-01  7.1635E-01
            -5.7780E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2640.33185211898        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  1.0000E+00  1.0804E+00  2.4863E+02  1.0379E+00  1.7524E+00  1.0352E+00  6.7498E-01  1.4778E-02  1.0224E+00  1.9513E+00
             2.8314E+00
 PARAMETER:  1.0004E-01  1.7731E-01  5.6160E+00  1.3715E-01  6.6100E-01  1.3455E-01 -2.9308E-01 -4.1146E+00  1.2213E-01  7.6852E-01
             1.1408E+00
 GRADIENT:   2.9253E+00 -1.6391E+01 -3.6649E-01 -2.0086E+01 -8.1489E+00 -3.4442E-01  2.6747E+00 -7.9732E-07  4.8826E+00 -4.0345E+00
             1.2120E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2644.13187658327        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      575
 NPARAMETR:  1.0303E+00  1.0974E+00  1.9821E+04  1.0538E+00  1.8887E+00  1.0663E+00  5.0425E-01  1.0000E-02  1.0532E+00  2.0471E+00
             2.8442E+00
 PARAMETER:  1.2982E-01  1.9294E-01  9.9945E+00  1.5245E-01  7.3591E-01  1.6422E-01 -5.8468E-01 -4.6864E+00  1.5179E-01  8.1644E-01
             1.1453E+00
 GRADIENT:  -3.1395E+00 -7.9129E-01 -6.6461E-03 -2.1451E+00 -1.0168E+00 -4.4966E-01 -6.2040E-01  0.0000E+00 -4.9095E-01 -4.4829E-01
             1.3430E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2644.21406089336        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      752            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0318E+00  1.0716E+00  1.4176E+04  1.0763E+00  1.8926E+00  1.0680E+00  6.1740E-01  1.0000E-02  1.0132E+00  2.0520E+00
             2.8424E+00
 PARAMETER:  1.3131E-01  1.6915E-01  9.6593E+00  1.7351E-01  7.3794E-01  1.6578E-01 -3.8224E-01 -4.5967E+00  1.1309E-01  8.1883E-01
             1.1446E+00
 GRADIENT:   8.2174E+01  1.3428E+01 -9.4877E-03  2.6294E+01  2.6489E+01  1.4639E+01  1.1359E+00  0.0000E+00  2.4510E+00  1.1023E+01
             2.0554E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2644.22060878452        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.0319E+00  1.0713E+00  6.8256E+04  1.0759E+00  1.8894E+00  1.0674E+00  6.2894E-01  1.0000E-02  1.0105E+00  2.0508E+00
             2.8424E+00
 PARAMETER:  1.3143E-01  1.6887E-01  1.1231E+01  1.7316E-01  7.3627E-01  1.6523E-01 -3.6372E-01 -4.5967E+00  1.1046E-01  8.1822E-01
             1.1446E+00
 GRADIENT:  -4.5575E+00 -1.0880E+00 -2.0625E-03 -1.1457E+00 -1.0906E+00 -4.2026E-01  1.2606E-01  0.0000E+00  9.0588E-01  1.1782E-01
             5.6576E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2644.22709438638        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1030
 NPARAMETR:  1.0323E+00  1.0730E+00  1.9436E+05  1.0769E+00  1.8927E+00  1.0681E+00  6.2956E-01  1.0000E-02  1.0069E+00  2.0518E+00
             2.8424E+00
 PARAMETER:  1.3180E-01  1.7044E-01  1.2277E+01  1.7412E-01  7.3799E-01  1.6584E-01 -3.6273E-01 -4.5967E+00  1.0692E-01  8.1872E-01
             1.1446E+00
 GRADIENT:  -7.3650E+01  1.6246E+01 -7.3057E-04  2.0233E+01 -2.7033E+00 -5.4603E+00 -1.4457E+00  0.0000E+00 -3.7944E+00  1.3018E+00
             1.8296E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2644.22838163280        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1223             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0318E+00  1.0725E+00  8.2520E+05  1.0768E+00  1.8928E+00  1.0680E+00  6.3341E-01  1.0004E-02  1.0029E+00  2.0517E+00
             2.8422E+00
 PARAMETER:  1.3129E-01  1.6995E-01  1.3723E+01  1.7403E-01  7.3808E-01  1.6583E-01 -3.5663E-01 -4.5048E+00  1.0285E-01  8.1869E-01
             1.1446E+00
 GRADIENT:   2.6402E+01  2.5802E+01 -1.5486E-04  4.2511E+01  2.4853E+01  1.4421E+01  4.7711E-01  1.3163E-05  2.3816E-01  1.1099E+01
             2.1279E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2644.22889461103        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  1.0318E+00  1.0727E+00  8.8009E+05  1.0772E+00  1.8927E+00  1.0680E+00  6.3681E-01  1.0948E-02  1.0034E+00  2.0517E+00
             2.8421E+00
 PARAMETER:  1.3131E-01  1.7018E-01  1.3788E+01  1.7432E-01  7.3802E-01  1.6582E-01 -3.5128E-01 -4.4146E+00  1.0335E-01  8.1867E-01
             1.1446E+00
 GRADIENT:  -5.8907E+01  9.8133E+00 -1.6212E-04  1.3441E+01 -1.9566E+00  2.3697E+00 -2.6284E-01 -3.2445E-06 -4.7259E-01  5.1335E-01
             1.5925E+00

0ITERATION NO.:   63    OBJECTIVE VALUE:  -2644.22905789630        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1511
 NPARAMETR:  1.0318E+00  1.0725E+00  9.4290E+05  1.0771E+00  1.8928E+00  1.0679E+00  6.4079E-01  1.3845E-02  1.0034E+00  2.0517E+00
             2.8422E+00
 PARAMETER:  1.3131E-01  1.7002E-01  1.3857E+01  1.7427E-01  7.3806E-01  1.6572E-01 -3.4505E-01 -4.1798E+00  1.0336E-01  8.1866E-01
             1.1446E+00
 GRADIENT:  -2.9768E-02 -1.2721E-01 -1.8015E-04 -2.0980E-01  2.4706E-03 -3.0136E-02  1.2509E-02  3.2747E-04  1.3861E-01 -9.5352E-03
             5.6369E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1511
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0996E-03 -2.0093E-02  9.5479E-10  1.5571E-04 -1.3566E-02
 SE:             2.9346E-02  1.1539E-02  6.0960E-10  2.5442E-02  2.7136E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7011E-01  8.1614E-02  1.1729E-01  9.9512E-01  6.1713E-01

 ETASHRINKSD(%)  1.6865E+00  6.1345E+01  1.0000E+02  1.4764E+01  9.0905E+00
 ETASHRINKVR(%)  3.3446E+00  8.5058E+01  1.0000E+02  2.7349E+01  1.7355E+01
 EBVSHRINKSD(%)  1.6288E+00  6.1671E+01  1.0000E+02  1.4513E+01  7.0660E+00
 EBVSHRINKVR(%)  3.2311E+00  8.5309E+01  1.0000E+02  2.6920E+01  1.3633E+01
 RELATIVEINF(%)  9.6678E+01  4.5792E-01  0.0000E+00  2.2775E+00  8.1778E+01
 EPSSHRINKSD(%)  1.5551E+01
 EPSSHRINKVR(%)  2.8684E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          884
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1624.6833267058612     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2644.2290578963043     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1019.5457311904431     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.35
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2644.229       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.07E+00  9.43E+05  1.08E+00  1.89E+00  1.07E+00  6.41E-01  1.38E-02  1.00E+00  2.05E+00  2.84E+00
 


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
+        9.02E+02
 
 TH 2
+       -9.35E+00  3.58E+02
 
 TH 3
+       -7.39E-08 -2.41E-08  1.21E-16
 
 TH 4
+       -1.49E+01  4.85E+02  4.38E-08  7.03E+02
 
 TH 5
+       -3.02E+00 -1.08E+01  1.42E-09 -1.51E+01  6.42E+01
 
 TH 6
+        1.59E+01 -4.99E+00 -9.00E-08 -7.27E+00 -5.18E-01  1.69E+02
 
 TH 7
+       -3.66E+00 -2.99E+01  6.18E-08 -1.92E+01  5.06E-01 -3.65E+00  1.82E+01
 
 TH 8
+       -2.67E+01  2.85E+01  3.93E-08  2.17E+01  1.22E+00 -1.37E+01  2.78E+01  5.12E+01
 
 TH 9
+        2.12E+01 -3.65E+01  1.10E-07  1.96E+01  1.77E+00 -8.30E-01  3.16E+01  2.71E+00  1.27E+02
 
 TH10
+        1.05E-01  1.20E+00  1.28E-08 -1.09E+00 -4.34E+00 -8.85E-01 -4.37E-01  1.58E-02 -2.19E+00  3.31E+01
 
 TH11
+       -1.37E+01 -1.65E+01  2.84E-10 -1.98E+01  1.35E+00  1.91E+00  2.68E+00 -4.50E-01  6.55E+00  4.35E+00  1.45E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       54.874
Stop Time:
Wed Sep 29 04:35:52 CDT 2021
