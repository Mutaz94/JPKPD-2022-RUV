Wed Sep 29 18:53:31 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat28.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1622.23061569314        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0659E+02 -8.7943E+01 -2.6883E+01 -1.0672E+02  6.9251E+01  3.4662E+01  4.7462E+00  7.2621E+00 -2.0786E+01  1.0875E+01
            -3.9143E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1636.68646736480        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9541E-01  1.0198E+00  1.0300E+00  1.0984E+00  9.6190E-01  9.8548E-01  9.5210E-01  9.7221E-01  1.0503E+00  9.3034E-01
             1.1207E+00
 PARAMETER:  9.5400E-02  1.1957E-01  1.2960E-01  1.9383E-01  6.1150E-02  8.5369E-02  5.0919E-02  7.1813E-02  1.4904E-01  2.7793E-02
             2.1399E-01
 GRADIENT:   2.3756E+01 -2.3037E+00 -3.7230E+00 -7.5787E+00  1.0426E+01 -7.2297E-01  8.1138E+00  3.5313E+00  8.3361E-01  6.5063E+00
             9.6689E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1638.31770844407        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.8965E-01  8.9589E-01  9.1142E-01  1.1928E+00  8.3894E-01  1.0046E+00  6.9847E-01  7.7325E-01  1.0384E+00  8.5836E-01
             1.1048E+00
 PARAMETER:  8.9598E-02 -9.9382E-03  7.2537E-03  2.7632E-01 -7.5614E-02  1.0460E-01 -2.5887E-01 -1.5715E-01  1.3772E-01 -5.2732E-02
             1.9962E-01
 GRADIENT:   8.8982E+00  1.9885E+01 -5.0503E+00  3.6739E+01 -6.0500E+00  6.6135E+00 -1.5635E+00  1.7438E+00  1.2156E+00  3.2522E+00
             3.7784E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1639.33332819910        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8443E-01  8.0114E-01  9.0526E-01  1.2281E+00  8.0108E-01  9.8307E-01  9.4756E-01  6.9439E-01  9.6238E-01  8.0743E-01
             1.0934E+00
 PARAMETER:  8.4306E-02 -1.2171E-01  4.6408E-04  3.0545E-01 -1.2179E-01  8.2921E-02  4.6140E-02 -2.6472E-01  6.1651E-02 -1.1390E-01
             1.8925E-01
 GRADIENT:  -1.1204E+00  3.4573E+00  5.0692E-02  4.9717E-01  1.4256E-01 -1.6747E+00  3.7651E-01 -1.1428E-01 -6.2842E-01 -3.4329E-01
            -1.0066E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1640.19083309696        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.8142E-01  4.7682E-01  1.0012E+00  1.4366E+00  7.4314E-01  9.8686E-01  9.9124E-01  7.6751E-01  8.6214E-01  8.0819E-01
             1.0956E+00
 PARAMETER:  8.1241E-02 -6.4061E-01  1.0124E-01  4.6226E-01 -1.9687E-01  8.6774E-02  9.1197E-02 -1.6460E-01 -4.8341E-02 -1.1296E-01
             1.9129E-01
 GRADIENT:   1.6886E+00  7.1287E+00  6.9834E+00  1.6828E+01 -1.1450E+01  1.1436E+00 -2.3620E-01 -3.6895E-01 -2.0004E-01  1.3458E-02
            -8.1054E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1640.43097584886        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.7778E-01  3.4663E-01  1.0382E+00  1.5130E+00  7.2972E-01  9.7975E-01  1.0098E+00  8.2619E-01  8.2106E-01  7.9505E-01
             1.0930E+00
 PARAMETER:  7.7534E-02 -9.5950E-01  1.3753E-01  5.1407E-01 -2.1509E-01  7.9541E-02  1.0979E-01 -9.0925E-02 -9.7155E-02 -1.2935E-01
             1.8894E-01
 GRADIENT:  -9.7967E-01  3.6016E+00  3.2152E+00  1.2733E+01 -3.0872E+00 -1.0059E+00 -2.7900E-01 -4.5059E-01 -1.2582E+00 -1.6290E+00
            -1.4656E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1640.47862829143        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.7656E-01  2.7471E-01  1.0393E+00  1.5557E+00  7.1148E-01  9.8079E-01  1.0258E+00  8.3656E-01  7.9963E-01  7.9604E-01
             1.0950E+00
 PARAMETER:  7.6281E-02 -1.1921E+00  1.3857E-01  5.4192E-01 -2.4041E-01  8.0602E-02  1.2544E-01 -7.8452E-02 -1.2360E-01 -1.2811E-01
             1.9079E-01
 GRADIENT:  -7.6153E-01  3.6319E+00  5.3895E+00  1.6172E+01 -8.3537E+00 -2.5779E-01 -2.3727E-01 -3.3983E-01 -1.6014E+00 -1.3502E-01
            -5.0688E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1640.59049236799        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1236
 NPARAMETR:  9.7489E-01  1.8159E-01  1.0290E+00  1.6038E+00  6.8599E-01  9.8143E-01  1.1026E+00  8.4675E-01  7.7507E-01  7.9231E-01
             1.0972E+00
 PARAMETER:  7.4569E-02 -1.6060E+00  1.2859E-01  5.7240E-01 -2.7690E-01  8.1255E-02  1.9763E-01 -6.6355E-02 -1.5480E-01 -1.3281E-01
             1.9278E-01
 GRADIENT:  -2.6300E-01  1.9047E+00  4.1793E+00  8.6963E+00 -8.6066E+00  4.3811E-01 -1.1604E-01  9.9267E-02 -7.2913E-01  1.4395E+00
             6.6645E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1640.63847656721        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1412
 NPARAMETR:  9.7365E-01  1.2558E-01  1.0366E+00  1.6352E+00  6.7788E-01  9.8055E-01  1.2126E+00  8.6966E-01  7.5933E-01  7.8433E-01
             1.0972E+00
 PARAMETER:  7.3294E-02 -1.9748E+00  1.3591E-01  5.9175E-01 -2.8879E-01  8.0360E-02  2.9277E-01 -3.9654E-02 -1.7531E-01 -1.4293E-01
             1.9280E-01
 GRADIENT:  -1.5250E-01  1.1362E+00  2.7254E+00  5.6250E+00 -6.0090E+00  3.8181E-01 -6.3889E-02  3.7814E-02 -6.3022E-01  1.0744E+00
             5.4747E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1640.70622547610        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1597             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7343E-01  1.0224E-01  1.0414E+00  1.6423E+00  6.7723E-01  9.7952E-01  2.0888E+00  8.8045E-01  7.5499E-01  7.7775E-01
             1.0963E+00
 PARAMETER:  7.3073E-02 -2.1805E+00  1.4055E-01  5.9612E-01 -2.8974E-01  7.9305E-02  8.3657E-01 -2.7324E-02 -1.8105E-01 -1.5135E-01
             1.9194E-01
 GRADIENT:   2.9373E+02  7.1494E+00  3.6010E+00  6.7589E+02  2.5354E+01  2.4689E+01  7.2255E-01  1.5593E-01  1.5549E+01  1.2164E+00
             1.7754E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1640.70962449119        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1775
 NPARAMETR:  9.7198E-01  9.7992E-02  1.0425E+00  1.6446E+00  6.7626E-01  9.7872E-01  1.9681E+00  8.8403E-01  7.5078E-01  7.7552E-01
             1.0959E+00
 PARAMETER:  7.1581E-02 -2.2229E+00  1.4167E-01  5.9750E-01 -2.9118E-01  7.8495E-02  7.7707E-01 -2.3265E-02 -1.8664E-01 -1.5422E-01
             1.9162E-01
 GRADIENT:  -2.3382E+00  1.5467E-01  2.4812E-01 -1.1009E+01  7.2473E-01 -1.8552E-01 -1.5059E-02 -1.2419E-01 -2.4749E-02  1.4454E-01
            -1.3732E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1640.71309474970        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1939             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7291E-01  9.7462E-02  1.0422E+00  1.6447E+00  6.7560E-01  9.7930E-01  2.0449E+00  8.8541E-01  7.4964E-01  7.7197E-01
             1.0959E+00
 PARAMETER:  7.2537E-02 -2.2283E+00  1.4132E-01  5.9757E-01 -2.9215E-01  7.9079E-02  8.1535E-01 -2.1706E-02 -1.8816E-01 -1.5881E-01
             1.9157E-01
 GRADIENT:   2.9292E+02  6.7150E+00  5.0102E+00  6.8002E+02  2.3759E+01  2.4681E+01  5.8474E-01  6.0541E-02  1.4048E+01  6.8372E-01
             1.3904E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1640.71494919996        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2125             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7335E-01  9.6500E-02  1.0417E+00  1.6450E+00  6.7519E-01  9.7928E-01  2.1314E+00  8.8612E-01  7.4995E-01  7.7233E-01
             1.0959E+00
 PARAMETER:  7.2984E-02 -2.2382E+00  1.4083E-01  5.9777E-01 -2.9275E-01  7.9058E-02  8.5679E-01 -2.0899E-02 -1.8775E-01 -1.5834E-01
             1.9162E-01
 GRADIENT:   2.9394E+02  6.7021E+00  4.8712E+00  6.8037E+02  2.3809E+01  2.4636E+01  6.6428E-01  1.8826E-01  1.4386E+01  8.5358E-01
             1.4766E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1640.71563317554        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2299
 NPARAMETR:  9.7316E-01  9.5215E-02  1.0412E+00  1.6456E+00  6.7486E-01  9.7919E-01  2.1948E+00  8.8669E-01  7.4994E-01  7.7252E-01
             1.0960E+00
 PARAMETER:  7.2790E-02 -2.2516E+00  1.4036E-01  5.9813E-01 -2.9325E-01  7.8968E-02  8.8608E-01 -2.0264E-02 -1.8777E-01 -1.5809E-01
             1.9171E-01
 GRADIENT:   6.0052E-01  1.3217E-01  2.7256E-01 -1.1991E+01  5.6217E-01  1.7199E-02  9.3188E-03  4.1492E-02  1.9774E-01  7.6258E-02
             2.9896E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1640.71838771730        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     2459
 NPARAMETR:  9.7291E-01  8.9497E-02  1.0365E+00  1.6496E+00  6.7174E-01  9.7913E-01  2.2972E+00  8.8471E-01  7.4718E-01  7.6852E-01
             1.0958E+00
 PARAMETER:  7.2539E-02 -2.3136E+00  1.3581E-01  6.0053E-01 -2.9788E-01  7.8910E-02  9.3169E-01 -2.2497E-02 -1.9145E-01 -1.6329E-01
             1.9152E-01
 GRADIENT:   2.8564E-01  2.2049E-01 -4.2902E-01 -9.0388E+00  1.1024E+00  1.2077E-02  7.2496E-04 -1.4444E-02 -4.4736E-01 -2.4536E-01
            -1.7811E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1640.71904975969        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2638
 NPARAMETR:  9.7261E-01  8.0710E-02  1.0317E+00  1.6550E+00  6.6771E-01  9.7904E-01  2.5648E+00  8.8231E-01  7.4462E-01  7.6534E-01
             1.0959E+00
 PARAMETER:  7.2227E-02 -2.4169E+00  1.3121E-01  6.0379E-01 -3.0391E-01  7.8813E-02  1.0419E+00 -2.5207E-02 -1.9488E-01 -1.6743E-01
             1.9159E-01
 GRADIENT:   2.9615E-03  2.8857E-01 -5.2540E-01 -6.6374E+00  7.4954E-01  6.2222E-03  1.3303E-02 -7.5786E-02 -5.2011E-01 -3.5072E-01
            -2.2664E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1640.72580846730        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2822
 NPARAMETR:  9.7272E-01  7.4796E-02  1.0308E+00  1.6554E+00  6.6637E-01  9.7905E-01  2.5957E+00  8.8272E-01  7.4470E-01  7.6664E-01
             1.0962E+00
 PARAMETER:  7.2342E-02 -2.4930E+00  1.3032E-01  6.0404E-01 -3.0591E-01  7.8823E-02  1.0539E+00 -2.4752E-02 -1.9477E-01 -1.6573E-01
             1.9188E-01
 GRADIENT:   6.4984E-01  7.5871E-02 -8.3492E-01 -1.2945E+01  1.5987E+00  5.3452E-02  1.4381E-02  4.3313E-02  1.3521E-01 -9.4814E-02
             1.9888E-02

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1640.72618352043        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     2959
 NPARAMETR:  9.7285E-01  7.3695E-02  1.0305E+00  1.6559E+00  6.6612E-01  9.7901E-01  2.5359E+00  8.8303E-01  7.4459E-01  7.6709E-01
             1.0962E+00
 PARAMETER:  7.2477E-02 -2.4990E+00  1.3133E-01  6.0431E-01 -3.0678E-01  7.8804E-02  1.0246E+00 -2.5157E-02 -1.9498E-01 -1.6475E-01
             1.9186E-01
 GRADIENT:   4.7411E-03  3.8763E-02  7.4384E-01 -5.6551E-02 -1.0830E+00  5.8240E-03 -8.7680E-04 -3.1840E-02 -2.2327E-02  5.3173E-02
            -1.6401E-03
 NUMSIGDIG:         4.7         2.4         1.9         4.6         2.7         3.8         2.2         2.1         3.5         2.6
                    4.7

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2959
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2178E-04 -2.2572E-03 -1.8962E-02 -4.6152E-03 -2.3634E-02
 SE:             2.9813E-02  3.1013E-03  1.7617E-02  2.8988E-02  2.0742E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9406E-01  4.6673E-01  2.8176E-01  8.7350E-01  2.5453E-01

 ETASHRINKSD(%)  1.2327E-01  8.9610E+01  4.0982E+01  2.8881E+00  3.0511E+01
 ETASHRINKVR(%)  2.4638E-01  9.8921E+01  6.5168E+01  5.6927E+00  5.1713E+01
 EBVSHRINKSD(%)  5.0184E-01  8.9872E+01  4.2419E+01  3.2155E+00  2.9329E+01
 EBVSHRINKVR(%)  1.0012E+00  9.8974E+01  6.6845E+01  6.3277E+00  5.0057E+01
 RELATIVEINF(%)  9.3841E+01  4.0224E-02  3.6622E+00  5.7105E+00  2.6711E+00
 EPSSHRINKSD(%)  4.4028E+01
 EPSSHRINKVR(%)  6.8671E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1640.7261835204322     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -905.57535695669401     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.35
 Elapsed covariance  time in seconds:     5.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1640.726       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  7.44E-02  1.03E+00  1.66E+00  6.66E-01  9.79E-01  2.52E+00  8.82E-01  7.45E-01  7.67E-01  1.10E+00
 


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
 
         2.96E-02  2.64E-01  2.35E-01  1.63E-01  1.03E-01  5.99E-02  5.31E+00  2.79E-01  8.79E-02  1.61E-01  8.05E-02
 


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
+        8.76E-04
 
 TH 2
+        1.81E-03  6.96E-02
 
 TH 3
+       -1.13E-03 -1.93E-02  5.52E-02
 
 TH 4
+       -1.08E-03 -4.14E-02  1.66E-02  2.65E-02
 
 TH 5
+       -2.05E-04  5.61E-03  2.04E-02 -1.04E-03  1.05E-02
 
 TH 6
+       -3.50E-04 -3.48E-04 -2.47E-03 -1.97E-04 -1.23E-03  3.59E-03
 
 TH 7
+       -3.08E-02 -1.33E+00  2.78E-01  7.75E-01 -1.46E-01 -1.21E-02  2.82E+01
 
 TH 8
+       -1.97E-03 -3.54E-02  5.06E-02  2.57E-02  1.49E-02 -2.89E-03  6.55E-01  7.80E-02
 
 TH 9
+        5.81E-04  1.83E-02 -6.35E-03 -1.10E-02  1.01E-03 -4.38E-04 -3.05E-01 -9.92E-03  7.72E-03
 
 TH10
+        7.93E-04  1.21E-02  1.01E-02 -6.27E-03  7.50E-03  5.01E-04 -3.09E-01 -6.63E-03  2.19E-03  2.59E-02
 
 TH11
+       -4.91E-05 -1.88E-03 -2.64E-03  1.32E-03 -1.55E-03 -8.87E-05  4.97E-02 -3.33E-03  1.15E-04 -2.65E-03  6.49E-03
 
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
+        2.96E-02
 
 TH 2
+        2.31E-01  2.64E-01
 
 TH 3
+       -1.62E-01 -3.12E-01  2.35E-01
 
 TH 4
+       -2.23E-01 -9.63E-01  4.33E-01  1.63E-01
 
 TH 5
+       -6.73E-02  2.07E-01  8.43E-01 -6.22E-02  1.03E-01
 
 TH 6
+       -1.97E-01 -2.20E-02 -1.75E-01 -2.02E-02 -1.99E-01  5.99E-02
 
 TH 7
+       -1.96E-01 -9.47E-01  2.23E-01  8.96E-01 -2.67E-01 -3.81E-02  5.31E+00
 
 TH 8
+       -2.38E-01 -4.81E-01  7.72E-01  5.65E-01  5.19E-01 -1.73E-01  4.41E-01  2.79E-01
 
 TH 9
+        2.23E-01  7.89E-01 -3.08E-01 -7.69E-01  1.12E-01 -8.32E-02 -6.53E-01 -4.04E-01  8.79E-02
 
 TH10
+        1.66E-01  2.84E-01  2.68E-01 -2.39E-01  4.54E-01  5.19E-02 -3.62E-01 -1.47E-01  1.55E-01  1.61E-01
 
 TH11
+       -2.06E-02 -8.84E-02 -1.39E-01  1.01E-01 -1.87E-01 -1.84E-02  1.16E-01 -1.48E-01  1.62E-02 -2.05E-01  8.05E-02
 
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
+        1.40E+03
 
 TH 2
+       -1.46E+02  6.97E+02
 
 TH 3
+       -1.19E+02  1.98E+02  4.87E+02
 
 TH 4
+       -9.46E+01  4.88E+02 -1.26E+01  7.58E+02
 
 TH 5
+        3.59E+02 -5.22E+02 -9.79E+02 -1.26E+02  2.36E+03
 
 TH 6
+        1.72E+02  4.01E+01 -2.04E+01  1.77E+01  1.03E+02  3.34E+02
 
 TH 7
+       -1.10E+00  1.28E+01  1.08E+00  2.20E+00 -6.08E-01  2.17E+00  5.14E-01
 
 TH 8
+        1.52E+01 -9.24E+00 -4.34E+01 -1.73E+01 -1.44E+01  3.28E+00 -5.26E-01  5.40E+01
 
 TH 9
+       -3.66E+01 -2.26E+02  3.71E+01  1.45E+00 -9.40E+01  1.51E+00 -6.68E+00  2.03E+00  4.52E+02
 
 TH10
+       -6.32E+01  3.85E+01  2.98E-02  1.53E+01 -1.06E+02 -2.06E+01  8.59E-01  3.05E+01 -1.16E+00  7.78E+01
 
 TH11
+        1.85E+01 -2.53E+01 -8.20E+00 -6.84E+01  6.53E-01  6.78E+00 -1.50E-01  2.41E+01 -2.96E+01  2.30E+01  1.81E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       40.729
Stop Time:
Wed Sep 29 18:54:14 CDT 2021
