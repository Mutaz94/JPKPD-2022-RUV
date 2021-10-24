Sun Oct 24 04:13:08 CDT 2021
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
$DATA ../../../../data/SD4/D/dat7.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1628.18067935122        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1133E+02  4.5773E+01 -1.8559E+01  1.2380E+02  5.2425E+01 -2.7223E+00  9.9255E+00  1.2443E+00  2.1715E+01 -4.3995E+00
            -7.6119E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1632.43441149071        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.8950E-01  9.9882E-01  1.0164E+00  9.5075E-01  9.4375E-01  1.2877E+00  9.6879E-01  9.9903E-01  9.4712E-01  1.0127E+00
             1.0274E+00
 PARAMETER:  8.9444E-02  9.8819E-02  1.1627E-01  4.9496E-02  4.2111E-02  3.5283E-01  6.8288E-02  9.9031E-02  4.5675E-02  1.1260E-01
             1.2707E-01
 GRADIENT:   1.7114E+01 -1.2100E+01  2.9147E+01 -5.9428E+01 -5.1814E+01  3.6885E+01  3.6115E+00 -1.1357E+00 -1.7481E+00  4.4026E+00
             4.0815E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1633.41356784432        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.8727E-01  1.0044E+00  9.3351E-01  9.4996E-01  9.2234E-01  1.2594E+00  8.3811E-01  7.6050E-01  9.6914E-01  9.4507E-01
             1.0566E+00
 PARAMETER:  8.7188E-02  1.0434E-01  3.1202E-02  4.8663E-02  1.9156E-02  3.3066E-01 -7.6609E-02 -1.7378E-01  6.8652E-02  4.3501E-02
             1.5507E-01
 GRADIENT:   1.3384E+01 -1.0279E+01  2.5862E+01 -4.2185E+01 -3.3420E+01  2.8896E+01 -5.0146E+00 -3.6006E+00 -2.9110E+00 -9.5152E+00
             1.0336E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1636.39092209904        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.7918E-01  1.1815E+00  7.4127E-01  8.5705E-01  9.2811E-01  1.1739E+00  9.0170E-01  6.1166E-01  1.0021E+00  9.5501E-01
             1.0239E+00
 PARAMETER:  7.8956E-02  2.6680E-01 -1.9939E-01 -5.4259E-02  2.5390E-02  2.6036E-01 -3.4787E-03 -3.9159E-01  1.0207E-01  5.3970E-02
             1.2366E-01
 GRADIENT:  -1.1668E+00  8.6281E+00  1.3937E+00  7.7578E+00 -4.7883E+00  1.2753E+00  8.3705E-01 -4.3653E-02 -7.6797E-01  3.0528E-01
             2.1383E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1636.58482261639        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.8006E-01  1.3695E+00  6.5918E-01  7.3291E-01  9.8545E-01  1.1725E+00  7.9876E-01  6.3094E-01  1.1333E+00  9.6950E-01
             1.0191E+00
 PARAMETER:  7.9855E-02  4.1444E-01 -3.1676E-01 -2.1074E-01  8.5339E-02  2.5914E-01 -1.2469E-01 -3.6054E-01  2.2511E-01  6.9021E-02
             1.1887E-01
 GRADIENT:  -3.8615E-01  7.1337E+00  3.3156E+00  2.4456E+00 -5.3689E+00  5.3692E-01 -2.8039E-01 -1.6269E-01 -4.5560E-01 -7.2772E-01
            -2.9277E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1636.61262614616        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.8029E-01  1.4453E+00  6.0682E-01  6.8220E-01  1.0015E+00  1.1730E+00  7.7259E-01  6.0181E-01  1.1881E+00  9.7283E-01
             1.0201E+00
 PARAMETER:  8.0089E-02  4.6830E-01 -3.9952E-01 -2.8243E-01  1.0155E-01  2.5958E-01 -1.5800E-01 -4.0781E-01  2.7238E-01  7.2454E-02
             1.1988E-01
 GRADIENT:  -3.7024E-01  5.5098E+00  1.9681E+00  2.3033E+00 -3.9888E+00  6.0438E-01 -2.5410E-01  2.2902E-02 -3.3578E-01 -1.1142E-01
             3.2398E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1636.63161039362        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1062             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8145E-01  1.4551E+00  5.9393E-01  6.7148E-01  1.0069E+00  1.1762E+00  7.6968E-01  5.8504E-01  1.2006E+00  9.7337E-01
             1.0196E+00
 PARAMETER:  8.1273E-02  4.7506E-01 -4.2100E-01 -2.9827E-01  1.0686E-01  2.6229E-01 -1.6178E-01 -4.3607E-01  2.8285E-01  7.3010E-02
             1.1936E-01
 GRADIENT:   4.4976E+02  4.0042E+02  3.5552E+00  1.1329E+02  1.1869E+01  1.9955E+02  8.6516E+00  2.9556E-01  1.4459E+01  6.5574E-01
             1.3766E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1636.64200748186        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1245             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8157E-01  1.4559E+00  5.9201E-01  6.7192E-01  1.0046E+00  1.1758E+00  7.7088E-01  5.3390E-01  1.1993E+00  9.7525E-01
             1.0190E+00
 PARAMETER:  8.1395E-02  4.7563E-01 -4.2423E-01 -2.9762E-01  1.0461E-01  2.6194E-01 -1.6023E-01 -5.2755E-01  2.8177E-01  7.4940E-02
             1.1882E-01
 GRADIENT:   4.5035E+02  4.0482E+02  5.1574E+00  1.1385E+02  9.1549E+00  1.9924E+02  8.5108E+00  1.4489E-01  1.4260E+01  6.3067E-01
             8.3126E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1636.64522645243        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1430             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8167E-01  1.4551E+00  5.9010E-01  6.7165E-01  1.0033E+00  1.1769E+00  7.7186E-01  5.2759E-01  1.1983E+00  9.7363E-01
             1.0190E+00
 PARAMETER:  8.1503E-02  4.7509E-01 -4.2747E-01 -2.9801E-01  1.0333E-01  2.6293E-01 -1.5895E-01 -5.3944E-01  2.8088E-01  7.3277E-02
             1.1887E-01
 GRADIENT:   4.5040E+02  4.0283E+02  5.0902E+00  1.1342E+02  9.3026E+00  2.0026E+02  8.5075E+00  1.5709E-01  1.4206E+01  6.1057E-01
             9.0052E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1636.64697163537        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1613
 NPARAMETR:  9.8168E-01  1.4546E+00  5.8789E-01  6.7183E-01  1.0023E+00  1.1770E+00  7.7235E-01  5.2088E-01  1.1974E+00  9.7209E-01
             1.0190E+00
 PARAMETER:  8.1511E-02  4.7471E-01 -4.3121E-01 -2.9775E-01  1.0225E-01  2.6293E-01 -1.5832E-01 -5.5223E-01  2.8014E-01  7.1697E-02
             1.1884E-01
 GRADIENT:   1.8564E+00 -3.8220E+00 -2.7587E-01  2.0798E-02  2.2632E+00  1.9479E+00 -1.1017E-01  2.2266E-02  1.7593E-01 -1.5135E-01
            -1.6296E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1636.64959784352        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1794
 NPARAMETR:  9.8169E-01  1.4545E+00  5.8646E-01  6.7185E-01  1.0009E+00  1.1770E+00  7.7298E-01  5.0845E-01  1.1964E+00  9.7105E-01
             1.0190E+00
 PARAMETER:  8.1520E-02  4.7469E-01 -4.3365E-01 -2.9772E-01  1.0090E-01  2.6294E-01 -1.5750E-01 -5.7638E-01  2.7928E-01  7.0625E-02
             1.1880E-01
 GRADIENT:   1.8595E+00 -3.1436E+00  7.2181E-02 -5.4286E-02  1.6053E+00  1.9492E+00 -1.2328E-01  1.9762E-03  1.1196E-01 -2.1906E-01
            -8.8755E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1636.65156151901        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1975
 NPARAMETR:  9.8170E-01  1.4540E+00  5.8464E-01  6.7173E-01  9.9976E-01  1.1770E+00  7.7364E-01  5.0208E-01  1.1954E+00  9.7050E-01
             1.0190E+00
 PARAMETER:  8.1535E-02  4.7435E-01 -4.3676E-01 -2.9790E-01  9.9763E-02  2.6296E-01 -1.5664E-01 -5.8899E-01  2.7851E-01  7.0059E-02
             1.1886E-01
 GRADIENT:   1.8700E+00 -3.8672E+00 -1.2592E-01 -3.1224E-01  1.7775E+00  1.9510E+00 -8.0134E-02  1.7023E-02  1.3762E-01 -1.0343E-01
            -2.6061E-03

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1636.65161840530        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     2034
 NPARAMETR:  9.8156E-01  1.4550E+00  5.8482E-01  6.7184E-01  9.9942E-01  1.1758E+00  7.7382E-01  4.9913E-01  1.1952E+00  9.7090E-01
             1.0190E+00
 PARAMETER:  8.1391E-02  4.7500E-01 -4.3645E-01 -2.9773E-01  9.9424E-02  2.6191E-01 -1.5641E-01 -5.9488E-01  2.7835E-01  7.0467E-02
             1.1886E-01
 GRADIENT:  -2.4386E-01  1.1641E+00  3.9162E-01  5.6198E-01  5.8964E-01 -4.1828E-01 -7.4247E-02 -1.6949E-03  8.7859E-03 -7.2426E-02
            -6.1264E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2034
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6706E-04 -2.7450E-02 -1.4836E-02  1.7571E-02 -3.1874E-02
 SE:             2.9917E-02  2.1745E-02  6.0478E-03  2.3795E-02  2.2899E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9288E-01  2.0681E-01  1.4165E-02  4.6026E-01  1.6395E-01

 ETASHRINKSD(%)  1.0000E-10  2.7153E+01  7.9739E+01  2.0282E+01  2.3285E+01
 ETASHRINKVR(%)  1.0000E-10  4.6933E+01  9.5895E+01  3.6451E+01  4.1147E+01
 EBVSHRINKSD(%)  3.2465E-01  2.6828E+01  8.2320E+01  2.1304E+01  2.1738E+01
 EBVSHRINKVR(%)  6.4824E-01  4.6459E+01  9.6874E+01  3.8070E+01  3.8751E+01
 RELATIVEINF(%)  9.9288E+01  2.3427E+00  2.6506E-01  2.9739E+00  9.6286E+00
 EPSSHRINKSD(%)  4.4715E+01
 EPSSHRINKVR(%)  6.9435E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1636.6516184053005     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -901.50079184156232     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1636.652       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.45E+00  5.85E-01  6.72E-01  9.99E-01  1.18E+00  7.74E-01  4.99E-01  1.20E+00  9.71E-01  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       60.129
Stop Time:
Sun Oct 24 04:13:20 CDT 2021
