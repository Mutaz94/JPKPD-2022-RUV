Sat Sep 18 12:44:14 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1638.44768509101        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.9123E+01 -7.4175E+01 -7.1901E+00 -1.0745E+02  6.0034E+00  3.3805E+01 -1.7017E+00  4.4724E+00 -3.8597E+00 -9.5365E+00
            -5.4129E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1648.41907459994        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9628E-01  1.0690E+00  1.0942E+00  1.0438E+00  1.0622E+00  8.8950E-01  9.8268E-01  9.5635E-01  9.7240E-01  1.0881E+00
             1.1106E+00
 PARAMETER:  9.6273E-02  1.6675E-01  1.9005E-01  1.4285E-01  1.6036E-01 -1.7093E-02  8.2529E-02  5.5368E-02  7.2010E-02  1.8446E-01
             2.0489E-01
 GRADIENT:   4.3410E+01  1.2242E+01  3.9104E+00  8.5878E+00 -7.6220E+00 -9.6706E+00  3.5113E-02 -2.2952E+00 -4.9752E+00 -5.8157E+00
            -6.3388E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1650.22542701246        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8943E-01  1.0099E+00  1.6628E+00  1.1034E+00  1.2850E+00  9.0454E-01  7.0720E-01  1.7057E+00  1.0591E+00  1.3812E+00
             1.1202E+00
 PARAMETER:  8.9379E-02  1.0985E-01  6.0853E-01  1.9844E-01  3.5080E-01 -3.2592E-04 -2.4645E-01  6.3399E-01  1.5743E-01  4.2298E-01
             2.1355E-01
 GRADIENT:   3.0450E+01  1.1619E+01 -5.5750E+00  2.8998E+01  2.9429E-01 -1.5653E+00  3.9760E-03  4.0992E+00  5.2666E+00  5.8404E+00
             1.9069E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1650.62150232869        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      294
 NPARAMETR:  9.8863E-01  9.8593E-01  1.6643E+00  1.1107E+00  1.2730E+00  9.1216E-01  8.1758E-01  1.6330E+00  1.0107E+00  1.3297E+00
             1.1192E+00
 PARAMETER:  8.8566E-02  8.5829E-02  6.0942E-01  2.0497E-01  3.4139E-01  8.0552E-03 -1.0140E-01  5.9044E-01  1.1066E-01  3.8497E-01
             2.1257E-01
 GRADIENT:  -3.3016E+00  1.8404E+00 -2.1449E+00  2.5028E+00  2.8421E+00 -8.6198E-01  6.0597E-01  1.0726E+00  9.4020E-01 -3.3676E-01
             5.3283E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1651.10317466284        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      470
 NPARAMETR:  9.8723E-01  7.4508E-01  1.8936E+00  1.2674E+00  1.2383E+00  9.1059E-01  7.8989E-01  1.6799E+00  9.2145E-01  1.3277E+00
             1.1162E+00
 PARAMETER:  8.7150E-02 -1.9427E-01  7.3850E-01  3.3697E-01  3.1373E-01  6.3423E-03 -1.3586E-01  6.1876E-01  1.8189E-02  3.8343E-01
             2.0994E-01
 GRADIENT:  -1.2848E+00  2.2833E+00  1.2946E+00  9.2557E-01 -1.5404E+00 -6.3563E-01 -1.1558E-01 -6.3421E-01 -4.1994E-01 -2.0699E-01
            -5.0716E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1651.55235807325        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  9.8453E-01  4.5560E-01  2.2207E+00  1.4584E+00  1.2209E+00  9.0770E-01  7.8326E-01  1.8912E+00  8.1981E-01  1.3303E+00
             1.1144E+00
 PARAMETER:  8.4410E-02 -6.8614E-01  8.9784E-01  4.7732E-01  2.9956E-01  3.1579E-03 -1.4429E-01  7.3722E-01 -9.8684E-02  3.8540E-01
             2.0834E-01
 GRADIENT:  -6.3172E-01  1.3761E+00  1.1134E+00  5.6476E-01 -1.9175E+00 -5.8874E-01 -1.5921E-01 -1.2189E-01 -1.4913E+00 -5.7071E-01
            -8.4680E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.00903843524        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  9.8246E-01  2.3371E-01  2.5662E+00  1.6151E+00  1.2302E+00  9.0905E-01  8.1968E-01  2.1677E+00  7.5006E-01  1.3525E+00
             1.1157E+00
 PARAMETER:  8.2303E-02 -1.3537E+00  1.0424E+00  5.7943E-01  3.0718E-01  4.6474E-03 -9.8844E-02  8.7365E-01 -1.8760E-01  4.0194E-01
             2.0946E-01
 GRADIENT:   7.2798E-02  2.5629E+00 -7.6613E-01  1.8379E+01 -1.5036E+00  1.0527E+00 -6.6886E-02  4.3431E-01 -7.2661E-01  3.6620E-01
            -7.3558E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1652.49991352475        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  9.8201E-01  1.0364E-01  2.9925E+00  1.7049E+00  1.2691E+00  9.0596E-01  9.5630E-01  2.4569E+00  7.0647E-01  1.3829E+00
             1.1175E+00
 PARAMETER:  8.1845E-02 -2.1669E+00  1.1961E+00  6.3349E-01  3.3831E-01  1.2369E-03  5.5317E-02  9.9890E-01 -2.4748E-01  4.2421E-01
             2.1112E-01
 GRADIENT:   1.9438E+00  1.0808E+00  2.5607E-01  1.5449E+01  5.9901E-01  4.8358E-01 -2.2361E-02 -8.6649E-01 -1.5652E+00  3.8076E-01
             6.3893E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1652.75868032339        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  9.8047E-01  4.4947E-02  3.1490E+00  1.7385E+00  1.2696E+00  9.0300E-01  1.1683E+00  2.5585E+00  6.9253E-01  1.3750E+00
             1.1140E+00
 PARAMETER:  8.0272E-02 -3.0023E+00  1.2471E+00  6.5302E-01  3.3870E-01 -2.0311E-03  2.5559E-01  1.0394E+00 -2.6741E-01  4.1844E-01
             2.0794E-01
 GRADIENT:  -4.4160E-01  1.9855E-01  1.4743E+00  1.2559E-01 -9.4717E-01 -4.6469E-01 -3.4795E-03 -8.4089E-01  2.3025E-01 -4.3860E-01
            -5.7045E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1652.76188971537        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1373
 NPARAMETR:  9.8050E-01  4.4213E-02  3.1502E+00  1.7390E+00  1.2696E+00  9.0324E-01  1.1735E+00  2.5595E+00  6.9219E-01  1.3751E+00
             1.1140E+00
 PARAMETER:  8.0305E-02 -3.0187E+00  1.2475E+00  6.5329E-01  3.3867E-01 -1.7716E-03  2.6003E-01  1.0398E+00 -2.6789E-01  4.1852E-01
             2.0797E-01
 GRADIENT:  -3.1404E-01  1.9500E-01  1.4636E+00  5.3995E-02 -9.6942E-01 -3.5997E-01 -3.4569E-03 -8.2884E-01  1.6795E-01 -4.1899E-01
            -5.5384E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1652.77220472283        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1510
 NPARAMETR:  9.8021E-01  3.8917E-02  3.1521E+00  1.7384E+00  1.2693E+00  9.0385E-01  1.1740E+00  2.5582E+00  6.9166E-01  1.3784E+00
             1.1156E+00
 PARAMETER:  8.0017E-02 -3.1463E+00  1.2481E+00  6.5296E-01  3.3850E-01 -1.0969E-03  2.6045E-01  1.0393E+00 -2.6866E-01  4.2092E-01
             2.0936E-01
 GRADIENT:   3.1600E+01  1.7514E-01  2.1904E+00  8.5937E+01  6.7708E-01  2.6517E+00 -4.2576E-04 -4.1988E-01  3.1954E+00  4.5866E-01
             4.5534E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1652.78486585585        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1689            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8050E-01  3.8089E-02  3.1422E+00  1.7395E+00  1.2696E+00  9.0389E-01  1.1744E+00  2.5629E+00  6.9024E-01  1.3785E+00
             1.1149E+00
 PARAMETER:  8.0307E-02 -3.1678E+00  1.2449E+00  6.5360E-01  3.3868E-01 -1.0432E-03  2.6078E-01  1.0411E+00 -2.7071E-01  4.2101E-01
             2.0874E-01
 GRADIENT:   3.2473E+01  1.9649E-01  1.5491E+00  8.7807E+01  1.1412E+00  2.6715E+00 -9.7494E-04  4.4272E-02  2.5928E+00  4.2710E-01
             1.6665E-01

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1652.78599710357        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     1784
 NPARAMETR:  9.8037E-01  3.7983E-02  3.1380E+00  1.7393E+00  1.2695E+00  9.0382E-01  1.1747E+00  2.5640E+00  6.9018E-01  1.3786E+00
             1.1149E+00
 PARAMETER:  8.0107E-02 -3.1702E+00  1.2434E+00  6.5356E-01  3.3864E-01 -1.1758E-03  2.6079E-01  1.0417E+00 -2.7084E-01  4.2099E-01
             2.0872E-01
 GRADIENT:  -4.2434E-01  8.0020E+03 -4.0792E+04  7.7610E+04  1.4984E+05 -4.9700E-02 -4.7064E-03  4.8641E+04 -2.9690E-02 -6.0271E+04
            -2.4310E+05
 NUMSIGDIG:         2.5         3.3         3.3         3.3         3.3         2.6         2.4         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1784
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5016E-04 -1.5270E-03 -4.3610E-02 -8.4534E-03 -6.4418E-02
 SE:             2.9785E-02  7.9682E-04  1.7992E-02  2.8954E-02  1.9676E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9598E-01  5.5321E-02  1.5354E-02  7.7031E-01  1.0607E-03

 ETASHRINKSD(%)  2.1665E-01  9.7331E+01  3.9726E+01  3.0020E+00  3.4083E+01
 ETASHRINKVR(%)  4.3283E-01  9.9929E+01  6.3671E+01  5.9139E+00  5.6550E+01
 EBVSHRINKSD(%)  6.4678E-01  9.7563E+01  4.4810E+01  3.3564E+00  2.8789E+01
 EBVSHRINKVR(%)  1.2894E+00  9.9941E+01  6.9541E+01  6.6001E+00  4.9290E+01
 RELATIVEINF(%)  9.6897E+01  3.3576E-03  1.5120E+01  5.6162E+00  1.9788E+01
 EPSSHRINKSD(%)  4.4189E+01
 EPSSHRINKVR(%)  6.8851E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.7859971035689     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.63517053983071     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.81
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.786       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  3.80E-02  3.14E+00  1.74E+00  1.27E+00  9.04E-01  1.17E+00  2.56E+00  6.90E-01  1.38E+00  1.11E+00
 


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
+        1.43E+03
 
 TH 2
+        4.05E+03  8.74E+08
 
 TH 3
+       -1.28E+02  8.41E+02  8.33E+05
 
 TH 4
+        4.16E+02 -2.42E+03  1.55E+03  9.81E+06
 
 TH 5
+        1.14E+03 -7.69E+03 -1.17E+03  3.92E+03  6.86E+07
 
 TH 6
+        3.03E+00  3.88E+03 -1.21E+02  4.08E+02  1.09E+03  2.60E+02
 
 TH 7
+       -1.08E+01 -4.11E-01 -3.37E-02  7.58E-01  3.70E+00  1.31E+01 -2.80E-01
 
 TH 8
+        1.85E+02 -1.22E+03  3.29E+03 -1.13E+04  1.67E+03  1.76E+02  1.00E-01  1.77E+06
 
 TH 9
+        1.07E+01  5.49E+03 -1.71E+02  6.06E+02  1.57E+03 -2.74E+00  2.54E+00  2.55E+02  3.69E+02
 
 TH10
+       -8.49E+02 -1.81E+08  2.07E+03 -7.10E+03 -1.88E+04 -8.09E+02 -1.04E+00 -3.02E+03 -1.16E+03  3.77E+07
 
 TH11
+       -2.13E+03  1.41E+04  1.57E+03 -5.40E+03 -1.42E+04 -2.01E+03 -2.04E+00 -2.29E+03 -2.88E+03  3.47E+04  2.34E+08
 
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
 #CPUT: Total CPU Time in Seconds,       29.247
Stop Time:
Sat Sep 18 12:44:45 CDT 2021
