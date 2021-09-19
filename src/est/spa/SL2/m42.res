Sat Sep 18 12:15:31 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1700.62526226081        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.1264E+01  2.2131E+01 -3.3029E+01  7.3949E+01  7.1358E+01  6.6426E-02  6.8529E+00  3.6225E+00  1.0203E+00 -1.4822E+00
            -1.6142E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.42193935428        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      176
 NPARAMETR:  1.0559E+00  9.5440E-01  1.0040E+00  9.8631E-01  9.2873E-01  9.9871E-01  9.4092E-01  9.7802E-01  9.9478E-01  9.5141E-01
             1.0678E+00
 PARAMETER:  1.5443E-01  5.3333E-02  1.0395E-01  8.6213E-02  2.6064E-02  9.8704E-02  3.9098E-02  7.7771E-02  9.4768E-02  5.0189E-02
             1.6556E-01
 GRADIENT:   8.0310E+00  2.4267E+00  6.5708E+00 -1.0651E+01 -6.8751E+00 -4.7645E-02  4.3697E+00  2.3962E+00  2.0719E-01  5.0850E+00
             1.0637E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1708.43337575715        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  1.0548E+00  8.3617E-01  9.3731E-01  1.0675E+00  8.5242E-01  1.0036E+00  8.2306E-01  8.2537E-01  9.7242E-01  8.7275E-01
             1.0239E+00
 PARAMETER:  1.5330E-01 -7.8927E-02  3.5255E-02  1.6530E-01 -5.9676E-02  1.0363E-01 -9.4724E-02 -9.1927E-02  7.2029E-02 -3.6107E-02
             1.2366E-01
 GRADIENT:   5.6608E+00  1.0348E+01  5.7938E-01  1.9351E+01  6.5539E-01  1.7094E+00 -2.0179E+00 -5.8561E-01 -6.7300E-01 -3.0078E+00
            -8.0279E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1708.80491558912        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  1.0536E+00  8.4446E-01  8.6327E-01  1.0484E+00  8.2084E-01  9.9951E-01  9.8819E-01  7.2419E-01  9.4372E-01  8.3986E-01
             1.0420E+00
 PARAMETER:  1.5221E-01 -6.9059E-02 -4.7026E-02  1.4722E-01 -9.7426E-02  9.9510E-02  8.8117E-02 -2.2270E-01  4.2077E-02 -7.4522E-02
             1.4110E-01
 GRADIENT:   1.5714E+00  1.1661E+00  2.3514E-01  1.5068E+00 -1.1516E-02 -6.6448E-02 -4.9700E-02 -6.1348E-02 -3.4177E-01 -2.2633E-01
            -1.2496E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1708.81872138908        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.0516E+00  7.6635E-01  9.2415E-01  1.0994E+00  8.2028E-01  9.9833E-01  1.0013E+00  7.7230E-01  9.2323E-01  8.5707E-01
             1.0427E+00
 PARAMETER:  1.5028E-01 -1.6612E-01  2.1123E-02  1.9474E-01 -9.8115E-02  9.8328E-02  1.0135E-01 -1.5838E-01  2.0126E-02 -5.4235E-02
             1.4180E-01
 GRADIENT:  -1.6141E-01  1.0053E+00  5.2552E-01  1.3608E+00 -1.0955E+00 -3.2059E-02  1.8650E-02  1.3417E-02  2.0495E-01 -5.3380E-02
             2.4697E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1708.83879455717        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0503E+00  6.3849E-01  1.0281E+00  1.1817E+00  8.2440E-01  9.9587E-01  9.8560E-01  8.5065E-01  8.9076E-01  8.8750E-01
             1.0416E+00
 PARAMETER:  1.4905E-01 -3.4864E-01  1.2776E-01  2.6698E-01 -9.3098E-02  9.5863E-02  8.5493E-02 -6.1755E-02 -1.5678E-02 -1.9347E-02
             1.4073E-01
 GRADIENT:   2.0139E+00  2.4813E-01  2.8642E-01  3.5697E-01  2.5158E-02 -1.0561E-01  1.3660E-01 -2.4021E-02  3.8624E-01 -1.9197E-01
            -1.9716E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1708.86520142352        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0466E+00  5.2579E-01  1.1047E+00  1.2530E+00  8.2180E-01  9.9446E-01  8.9186E-01  9.0578E-01  8.6275E-01  9.1066E-01
             1.0419E+00
 PARAMETER:  1.4551E-01 -5.4285E-01  1.9959E-01  3.2557E-01 -9.6263E-02  9.4441E-02 -1.4451E-02  1.0428E-03 -4.7630E-02  6.4183E-03
             1.4107E-01
 GRADIENT:  -1.3976E+00  6.2817E-01  6.0334E-01  1.2404E+00 -1.3734E+00  1.0195E-01  2.3373E-02 -8.4060E-02 -2.9828E-02  2.7854E-01
            -1.4902E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1708.87094309854        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1236
 NPARAMETR:  1.0464E+00  4.7846E-01  1.1379E+00  1.2821E+00  8.2200E-01  9.9336E-01  8.1662E-01  9.3800E-01  8.5156E-01  9.1546E-01
             1.0419E+00
 PARAMETER:  1.4531E-01 -6.3718E-01  2.2919E-01  3.4853E-01 -9.6020E-02  9.3342E-02 -1.0258E-01  3.5997E-02 -6.0691E-02  1.1676E-02
             1.4108E-01
 GRADIENT:   1.5024E-01  8.2434E-04 -1.6034E-02  4.6513E-02  8.9193E-02 -7.0248E-03  1.0633E-02  5.7844E-03 -1.3025E-02 -3.5730E-02
            -9.7034E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1708.87094532200        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1412
 NPARAMETR:  1.0463E+00  4.7700E-01  1.1386E+00  1.2830E+00  8.2185E-01  9.9335E-01  8.1366E-01  9.3839E-01  8.5125E-01  9.1574E-01
             1.0419E+00
 PARAMETER:  1.4525E-01 -6.4023E-01  2.2978E-01  3.4921E-01 -9.6193E-02  9.3328E-02 -1.0621E-01  3.6407E-02 -6.1049E-02  1.1976E-02
             1.4108E-01
 GRADIENT:   8.3114E-02 -1.3889E-02 -1.8946E-02 -2.6770E-03  6.8362E-02 -2.9716E-03  1.0289E-02  4.6573E-03 -1.9437E-03 -1.7147E-02
            -4.8181E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1708.87094655858        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1588
 NPARAMETR:  1.0463E+00  4.7611E-01  1.1390E+00  1.2835E+00  8.2178E-01  9.9334E-01  8.1100E-01  9.3866E-01  8.5110E-01  9.1595E-01
             1.0419E+00
 PARAMETER:  1.4522E-01 -6.4210E-01  2.3019E-01  3.4962E-01 -9.6280E-02  9.3317E-02 -1.0949E-01  3.6698E-02 -6.1220E-02  1.2202E-02
             1.4109E-01
 GRADIENT:   3.7893E-02 -1.9276E-02 -1.8124E-02 -2.4699E-02  4.6796E-02 -1.7960E-03  9.2602E-03  3.8415E-03  4.3636E-03 -5.1469E-03
            -1.4760E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1708.87094842040        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1764
 NPARAMETR:  1.0462E+00  4.7506E-01  1.1397E+00  1.2842E+00  8.2175E-01  9.9333E-01  8.0634E-01  9.3909E-01  8.5100E-01  9.1627E-01
             1.0419E+00
 PARAMETER:  1.4517E-01 -6.4431E-01  2.3077E-01  3.5011E-01 -9.6318E-02  9.3311E-02 -1.1525E-01  3.7161E-02 -6.1340E-02  1.2559E-02
             1.4109E-01
 GRADIENT:  -1.3393E-02 -2.3604E-02 -1.4668E-02 -4.9583E-02  1.6569E-02  2.1934E-03  6.8784E-03  2.4555E-03  1.1162E-02  9.1743E-03
             2.5206E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1708.87095861110        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1944
 NPARAMETR:  1.0462E+00  4.7445E-01  1.1407E+00  1.2845E+00  8.2197E-01  9.9333E-01  7.9702E-01  9.3988E-01  8.5125E-01  9.1682E-01
             1.0419E+00
 PARAMETER:  1.4514E-01 -6.4560E-01  2.3163E-01  3.5041E-01 -9.6048E-02  9.3304E-02 -1.2688E-01  3.8001E-02 -6.1047E-02  1.3156E-02
             1.4109E-01
 GRADIENT:  -5.8626E-02 -1.0910E-02 -2.6053E-03 -3.5205E-02 -2.1852E-02  3.3726E-03  2.1759E-04  1.5244E-04  9.5599E-03  1.6382E-02
             4.6198E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1708.87102026380        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2130             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0463E+00  4.7697E-01  1.1398E+00  1.2830E+00  8.2235E-01  9.9336E-01  7.9934E-01  9.3933E-01  8.5197E-01  9.1668E-01
             1.0419E+00
 PARAMETER:  1.4523E-01 -6.4030E-01  2.3084E-01  3.4920E-01 -9.5584E-02  9.3338E-02 -1.2397E-01  3.7412E-02 -6.0204E-02  1.3001E-02
             1.4108E-01
 GRADIENT:   6.6849E+01  7.1407E+00  5.3287E-01  4.9243E+01  8.5437E-01  4.2469E+00  6.9905E-02  2.2428E-02  1.1425E+00  6.8866E-02
             1.1440E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1708.87102083914        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2312
 NPARAMETR:  1.0463E+00  4.7700E-01  1.1398E+00  1.2830E+00  8.2235E-01  9.9336E-01  7.9971E-01  9.3927E-01  8.5196E-01  9.1667E-01
             1.0419E+00
 PARAMETER:  1.4522E-01 -6.4024E-01  2.3082E-01  3.4922E-01 -9.5593E-02  9.3339E-02 -1.2351E-01  3.7346E-02 -6.0220E-02  1.2988E-02
             1.4108E-01
 GRADIENT:   3.0219E-02  4.6456E-03  7.0764E-03  4.2675E-03  6.8097E-03  1.7710E-03 -2.2527E-04 -3.3485E-04  3.1760E-03 -3.4655E-04
            -9.4672E-04

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1708.87102083914        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     2336
 NPARAMETR:  1.0463E+00  4.7700E-01  1.1398E+00  1.2830E+00  8.2235E-01  9.9336E-01  7.9971E-01  9.3927E-01  8.5196E-01  9.1667E-01
             1.0419E+00
 PARAMETER:  1.4522E-01 -6.4024E-01  2.3082E-01  3.4922E-01 -9.5593E-02  9.3339E-02 -1.2351E-01  3.7346E-02 -6.0220E-02  1.2988E-02
             1.4108E-01
 GRADIENT:  -8.7863E-03  4.8587E-03  7.0919E-03  5.2952E-02  7.0091E-03  2.0920E-04  1.5942E-04  4.2770E-05  2.5932E-03  6.7075E-04
            -1.0366E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2336
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5790E-04 -1.1895E-02 -2.6023E-02 -1.7448E-03 -2.8374E-02
 SE:             2.9848E-02  7.1512E-03  1.5559E-02  2.8453E-02  2.2456E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9311E-01  9.6253E-02  9.4412E-02  9.5110E-01  2.0640E-01

 ETASHRINKSD(%)  5.6342E-03  7.6042E+01  4.7877E+01  4.6788E+00  2.4769E+01
 ETASHRINKVR(%)  1.1268E-02  9.4260E+01  7.2832E+01  9.1387E+00  4.3404E+01
 EBVSHRINKSD(%)  4.5683E-01  7.6815E+01  5.1364E+01  4.7127E+00  2.2410E+01
 EBVSHRINKVR(%)  9.1157E-01  9.4625E+01  7.6346E+01  9.2033E+00  3.9798E+01
 RELATIVEINF(%)  9.7229E+01  2.1943E-01  4.2787E+00  4.8601E+00  6.6239E+00
 EPSSHRINKSD(%)  4.4679E+01
 EPSSHRINKVR(%)  6.9396E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1708.8710208391442     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -973.72019427540602     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.34
 Elapsed covariance  time in seconds:     5.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1708.871       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  4.77E-01  1.14E+00  1.28E+00  8.22E-01  9.93E-01  8.00E-01  9.39E-01  8.52E-01  9.17E-01  1.04E+00
 


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
 
         3.09E-02  5.24E-01  3.02E-01  3.33E-01  8.59E-02  6.33E-02  3.88E-01  2.94E-01  1.78E-01  2.01E-01  6.84E-02
 


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
+        9.53E-04
 
 TH 2
+        1.74E-03  2.74E-01
 
 TH 3
+       -8.59E-04 -1.36E-01  9.09E-02
 
 TH 4
+       -1.21E-03 -1.74E-01  8.81E-02  1.11E-01
 
 TH 5
+        1.37E-04  2.18E-02 -5.58E-04 -1.30E-02  7.37E-03
 
 TH 6
+       -4.59E-05  4.68E-03 -4.51E-03 -2.83E-03 -3.21E-04  4.00E-03
 
 TH 7
+       -5.57E-05  7.01E-02 -7.05E-02 -4.78E-02 -1.49E-02  3.03E-03  1.51E-01
 
 TH 8
+       -7.15E-04 -1.05E-01  5.09E-02  6.52E-02 -1.13E-02 -5.27E-03  1.60E-02  8.63E-02
 
 TH 9
+        7.98E-04  8.58E-02 -4.34E-02 -5.46E-02  6.68E-03  1.12E-03  1.90E-02 -3.71E-02  3.16E-02
 
 TH10
+       -4.82E-04 -3.51E-02  3.67E-02  2.37E-02  7.97E-03 -1.50E-03 -7.19E-02 -4.26E-03 -1.06E-02  4.03E-02
 
 TH11
+       -6.38E-05  7.05E-03 -4.10E-03 -4.61E-03 -6.71E-05  4.95E-04  3.43E-03 -2.40E-03  1.51E-03 -2.16E-03  4.68E-03
 
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
+        3.09E-02
 
 TH 2
+        1.08E-01  5.24E-01
 
 TH 3
+       -9.23E-02 -8.63E-01  3.02E-01
 
 TH 4
+       -1.18E-01 -9.94E-01  8.76E-01  3.33E-01
 
 TH 5
+        5.16E-02  4.85E-01 -2.16E-02 -4.53E-01  8.59E-02
 
 TH 6
+       -2.35E-02  1.41E-01 -2.36E-01 -1.34E-01 -5.91E-02  6.33E-02
 
 TH 7
+       -4.64E-03  3.45E-01 -6.02E-01 -3.69E-01 -4.46E-01  1.23E-01  3.88E-01
 
 TH 8
+       -7.88E-02 -6.79E-01  5.74E-01  6.65E-01 -4.49E-01 -2.83E-01  1.40E-01  2.94E-01
 
 TH 9
+        1.45E-01  9.21E-01 -8.10E-01 -9.20E-01  4.38E-01  9.95E-02  2.75E-01 -7.11E-01  1.78E-01
 
 TH10
+       -7.77E-02 -3.34E-01  6.07E-01  3.54E-01  4.62E-01 -1.18E-01 -9.23E-01 -7.23E-02 -2.96E-01  2.01E-01
 
 TH11
+       -3.02E-02  1.97E-01 -1.99E-01 -2.02E-01 -1.14E-02  1.14E-01  1.29E-01 -1.19E-01  1.24E-01 -1.57E-01  6.84E-02
 
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
+        1.15E+03
 
 TH 2
+        6.43E+01  4.56E+02
 
 TH 3
+       -2.13E+01  1.09E+02  2.76E+02
 
 TH 4
+        1.28E+02  5.51E+02 -8.29E+01  9.52E+02
 
 TH 5
+        3.20E+01 -3.64E+02 -4.69E+02 -8.46E-01  1.32E+03
 
 TH 6
+       -1.24E+01 -9.56E+00  4.65E+01 -5.97E+01 -6.99E+01  3.18E+02
 
 TH 7
+        5.94E+01 -6.54E+00 -7.93E-01  2.88E+01  6.06E+01 -2.40E+01  6.79E+01
 
 TH 8
+       -1.47E+01  1.02E+01 -2.92E+01  1.68E+01  1.25E+00  3.87E+01 -2.52E+01  4.81E+01
 
 TH 9
+       -3.88E+01 -3.59E+01  2.74E-01  5.43E+01  1.15E+01  4.56E+01  1.11E+01  2.96E+01  2.55E+02
 
 TH10
+        1.03E+02  2.45E+01 -2.06E+01  6.39E+01 -3.59E+01 -1.81E+01  8.82E+01 -5.26E+00  2.76E+01  2.00E+02
 
 TH11
+        3.70E+01 -1.87E+01 -4.00E+01  4.34E+01  9.98E+01 -2.42E+01  1.61E+01  2.47E+00  4.01E+01  2.71E+01  2.43E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.227
Stop Time:
Sat Sep 18 12:16:08 CDT 2021
