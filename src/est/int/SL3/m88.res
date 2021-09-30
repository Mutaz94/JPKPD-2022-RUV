Wed Sep 29 04:51:28 CDT 2021
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
$DATA ../../../../data/int/SL3/dat88.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      978
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

 TOT. NO. OF OBS RECS:      878
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   442.587503238975        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3065E+02 -2.0477E+01  4.1278E+01  2.0150E+02  2.4816E+02  4.5560E+01 -1.4658E+02 -2.3637E+02 -1.7499E+02 -3.0014E+01
            -7.9129E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2323.83009722504        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0332E+00  1.4747E+00  8.9028E-01  8.8946E-01  1.0595E+00  8.3283E-01  1.1934E+00  1.1312E+00  9.3005E-01  9.3423E-01
             5.2213E+00
 PARAMETER:  1.3269E-01  4.8847E-01 -1.6217E-02 -1.7137E-02  1.5783E-01 -8.2928E-02  2.7680E-01  2.2329E-01  2.7481E-02  3.1965E-02
             1.7528E+00
 GRADIENT:  -8.2288E+01  6.7112E+01 -8.4991E+00  4.4306E+01 -8.3074E+00 -4.8162E+01  2.9959E+01  4.3849E+00  1.8270E+01  7.7370E+00
             7.7853E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2373.60750504202        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0124E+00  1.6014E+00  9.6418E-01  7.9940E-01  1.1337E+00  9.4570E-01  1.1856E+00  5.6135E+00  7.7587E-01  9.7801E-01
             4.6531E+00
 PARAMETER:  1.1230E-01  5.7090E-01  6.3526E-02 -1.2390E-01  2.2550E-01  4.4168E-02  2.7028E-01  1.8252E+00 -1.5376E-01  7.7762E-02
             1.6375E+00
 GRADIENT:  -7.2330E+01  9.4302E+01 -2.1438E+01  5.7877E+01 -9.4360E+01 -4.5469E+00  3.4177E+01  6.8258E+01  1.3266E+01 -1.5149E-01
             6.7660E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2536.77009561992        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.9255E-01  1.9139E+00  1.6338E+00  6.0463E-01  1.5495E+00  9.4624E-01  7.6898E-01  7.6970E+00  1.2295E+00  2.0620E+00
             2.9886E+00
 PARAMETER:  9.2524E-02  7.4915E-01  5.9089E-01 -4.0315E-01  5.3790E-01  4.4738E-02 -1.6269E-01  2.1408E+00  3.0664E-01  8.2370E-01
             1.1948E+00
 GRADIENT:  -2.8347E+01  2.0540E+02 -4.0816E+00  9.7382E+01 -9.0949E+00 -9.9906E+00  7.5051E-01  3.2562E+01  1.0363E+01  6.2975E+01
             6.5688E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2592.58402115840        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.9946E-01  1.3878E+00  2.2033E+00  8.1217E-01  1.3270E+00  9.6837E-01  1.1290E+00  3.1687E+00  6.3459E-01  1.2422E+00
             2.9104E+00
 PARAMETER:  9.9458E-02  4.2769E-01  8.8995E-01 -1.0804E-01  3.8293E-01  6.7862E-02  2.2133E-01  1.2533E+00 -3.5478E-01  3.1686E-01
             1.1683E+00
 GRADIENT:   1.8884E+00  3.3699E+00 -4.2659E+00 -3.5935E+01 -2.2284E+01  2.7268E+00  5.6228E+00 -6.7221E+00  2.4853E+00 -2.4805E+00
            -1.6058E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2596.65767418711        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      396
 NPARAMETR:  9.9886E-01  1.2539E+00  3.1919E+00  9.4215E-01  1.3662E+00  9.6175E-01  1.2721E+00  3.5131E+00  4.9785E-01  1.2350E+00
             2.9097E+00
 PARAMETER:  9.8855E-02  3.2624E-01  1.2606E+00  4.0406E-02  4.1202E-01  6.1002E-02  3.4069E-01  1.3565E+00 -5.9746E-01  3.1106E-01
             1.1681E+00
 GRADIENT:  -5.6809E+01 -2.2877E+00 -2.9445E+00  7.9570E+00 -1.7623E+01 -6.5391E+00  2.8800E-01 -1.1940E+01  4.8113E-01 -3.5679E+00
            -1.3046E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2597.43697187390        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      542
 NPARAMETR:  1.0240E+00  1.2510E+00  3.1691E+00  9.3327E-01  1.3655E+00  9.7572E-01  1.2664E+00  3.5282E+00  4.4092E-01  1.2352E+00
             2.9287E+00
 PARAMETER:  1.2375E-01  3.2395E-01  1.2534E+00  3.0937E-02  4.1155E-01  7.5419E-02  3.3619E-01  1.3608E+00 -7.1889E-01  3.1124E-01
             1.1746E+00
 GRADIENT:   6.6449E+01  1.7032E+01  1.6666E+00  6.3183E+00 -3.7127E+00  5.9063E+00  1.4769E+00 -6.7330E+00  6.6022E-01 -1.9707E+00
             1.9221E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2597.74693367872        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.0235E+00  1.3143E+00  3.1719E+00  9.0357E-01  1.3666E+00  9.7728E-01  1.2345E+00  3.5428E+00  4.8559E-01  1.2354E+00
             2.9289E+00
 PARAMETER:  1.2320E-01  3.7332E-01  1.2543E+00 -1.4007E-03  4.1236E-01  7.7015E-02  3.1065E-01  1.3649E+00 -6.2238E-01  3.1141E-01
             1.1746E+00
 GRADIENT:  -7.4327E-01  1.5647E+00  1.1618E+00 -9.0626E-01 -3.2607E+01  1.2337E+00  1.2554E+00 -1.2685E+01  1.0278E-01 -6.2109E+00
            -2.2210E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2598.43816503194        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0243E+00  1.3582E+00  3.1767E+00  8.8052E-01  1.3840E+00  9.7227E-01  1.1364E+00  3.7380E+00  5.7393E-01  1.2379E+00
             2.9382E+00
 PARAMETER:  1.2397E-01  4.0618E-01  1.2559E+00 -2.7239E-02  4.2495E-01  7.1880E-02  2.2787E-01  1.4185E+00 -4.5524E-01  3.1344E-01
             1.1778E+00
 GRADIENT:   1.4235E-01  1.0663E+00  4.1675E-01  5.2818E+00 -3.0453E+01 -7.5559E-01 -6.4000E+00 -7.6184E+00  6.0971E-01 -7.6704E+00
             6.3694E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2598.61142347544        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  1.0242E+00  1.3610E+00  3.1663E+00  8.7669E-01  1.3835E+00  9.7575E-01  1.1815E+00  3.7360E+00  5.0576E-01  1.2371E+00
             2.9460E+00
 PARAMETER:  1.2391E-01  4.0825E-01  1.2525E+00 -3.1601E-02  4.2465E-01  7.5450E-02  2.6682E-01  1.4180E+00 -5.8168E-01  3.1274E-01
             1.1805E+00
 GRADIENT:  -1.0615E-01  3.7286E+00  7.5011E-01  2.1854E+00 -3.2419E+01  6.2583E-01 -1.5707E+00 -8.4489E+00  1.4900E-01 -8.1488E+00
             1.1760E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2598.78102243108        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1237
 NPARAMETR:  1.0243E+00  1.3620E+00  3.1543E+00  8.7610E-01  1.3877E+00  9.7806E-01  1.1864E+00  3.7674E+00  4.5955E-01  1.2368E+00
             2.9291E+00
 PARAMETER:  1.2397E-01  4.0894E-01  1.2488E+00 -3.2280E-02  4.2768E-01  7.7811E-02  2.7093E-01  1.4264E+00 -6.7751E-01  3.1255E-01
             1.1747E+00
 GRADIENT:   3.6105E-01  5.8349E+00 -1.0317E-01  6.6378E+00 -3.0370E+01  1.3615E+00 -3.7532E+00 -8.3441E+00 -7.9642E-01 -8.7512E+00
            -1.6256E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2598.83078357766        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  1.0264E+00  1.3532E+00  3.1485E+00  8.7921E-01  1.3880E+00  9.7340E-01  1.2052E+00  3.7684E+00  4.7008E-01  1.2364E+00
             2.9335E+00
 PARAMETER:  1.2606E-01  4.0250E-01  1.2469E+00 -2.8729E-02  4.2783E-01  7.3041E-02  2.8667E-01  1.4266E+00 -6.5485E-01  3.1221E-01
             1.1762E+00
 GRADIENT:   5.1285E+00  1.9764E+00 -6.1468E-01  1.6024E+00 -2.7969E+01 -3.3738E-01 -4.0455E-02 -7.8633E+00 -7.2155E-02 -8.3684E+00
             2.6586E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2599.10758282294        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1603             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0203E+00  1.3584E+00  3.1190E+00  8.7472E-01  1.3962E+00  9.7766E-01  1.2055E+00  3.8255E+00  4.6618E-01  1.2512E+00
             2.9490E+00
 PARAMETER:  1.2014E-01  4.0630E-01  1.2375E+00 -3.3851E-02  4.3376E-01  7.7406E-02  2.8693E-01  1.4417E+00 -6.6318E-01  3.2412E-01
             1.1815E+00
 GRADIENT:   5.4210E+01  4.3823E+01  1.4426E+00  1.1305E+01 -7.0961E+00  6.3814E+00  4.9319E+00 -6.9276E-01  1.3404E+00 -4.5865E+00
             3.6706E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2599.59133334072        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1743
 NPARAMETR:  1.0215E+00  1.3539E+00  3.1331E+00  8.7114E-01  1.4615E+00  9.6878E-01  1.1978E+00  3.7972E+00  4.3973E-01  1.2491E+00
             2.9667E+00
 PARAMETER:  1.2130E-01  4.0299E-01  1.2420E+00 -3.7949E-02  4.7949E-01  6.8280E-02  2.8048E-01  1.4343E+00 -7.2160E-01  3.2246E-01
             1.1874E+00
 GRADIENT:  -7.1028E+00 -1.6051E+01 -8.6682E+00  1.7420E+00  1.0107E+01 -1.9141E+00 -1.4435E+00 -6.4329E+00 -2.0792E-02 -1.1800E+01
             2.7676E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2599.95919161748        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1916
 NPARAMETR:  1.0304E+00  1.3987E+00  3.1309E+00  8.5126E-01  1.4519E+00  9.7599E-01  1.1867E+00  3.7962E+00  4.5133E-01  1.2491E+00
             2.9664E+00
 PARAMETER:  1.2985E-01  4.3579E-01  1.2420E+00 -6.0033E-02  4.7193E-01  7.6218E-02  2.7317E-01  1.4347E+00 -6.9326E-01  3.2261E-01
             1.1868E+00
 GRADIENT:  -1.6512E+03  4.9734E+02  1.6592E+02  5.9757E+00 -4.1642E+00  9.1400E-01  2.3550E+00  1.4584E+02  9.2854E-02  6.5511E+02
            -1.6430E+02
 NUMSIGDIG:         2.3         2.3         2.3         1.0         1.7         1.3         1.1         2.3         1.5         2.3
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1916
 NO. OF SIG. DIGITS IN FINAL EST.:  1.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.8517E-03 -3.4950E-03 -2.8970E-02 -9.0392E-03 -2.6504E-02
 SE:             2.9206E-02  2.6144E-02  1.8371E-02  1.0707E-02  2.3108E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2222E-01  8.9365E-01  1.1481E-01  3.9856E-01  2.5140E-01

 ETASHRINKSD(%)  2.1556E+00  1.2416E+01  3.8455E+01  6.4129E+01  2.2584E+01
 ETASHRINKVR(%)  4.2647E+00  2.3290E+01  6.2123E+01  8.7132E+01  4.0067E+01
 EBVSHRINKSD(%)  2.0272E+00  1.1794E+01  4.4730E+01  6.6216E+01  2.3602E+01
 EBVSHRINKVR(%)  4.0134E+00  2.2197E+01  6.9453E+01  8.8586E+01  4.1634E+01
 RELATIVEINF(%)  9.5884E+01  9.7655E+00  1.1764E+01  1.3153E+00  2.8983E+01
 EPSSHRINKSD(%)  1.7018E+01
 EPSSHRINKVR(%)  3.1141E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          878
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1613.6560643074051     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2599.9591916174777     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -986.30312731007257     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.08
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2599.959       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.40E+00  3.13E+00  8.52E-01  1.45E+00  9.76E-01  1.19E+00  3.80E+00  4.52E-01  1.25E+00  2.96E+00
 


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
+        3.03E+05
 
 TH 2
+       -1.10E+01  1.49E+04
 
 TH 3
+        4.79E-01  9.91E+00  3.79E+02
 
 TH 4
+       -1.55E+03  7.70E+02  3.42E+01  8.59E+02
 
 TH 5
+        8.28E+02 -2.25E+02  1.26E+02  9.31E+04  1.19E+04
 
 TH 6
+        9.42E+00 -8.00E+00 -1.06E+00 -8.24E+00 -8.13E+04  1.92E+02
 
 TH 7
+       -5.26E+01  2.72E+01  1.92E+00 -6.81E+01 -2.45E+04 -2.26E-01  8.59E+01
 
 TH 8
+       -2.79E+02  9.99E+03 -1.55E+03 -1.18E+04 -1.04E+04 -1.06E-01  3.11E+03  1.28E+03
 
 TH 9
+       -1.29E+02  1.35E+01  4.40E+00 -4.17E+01 -2.53E+04 -3.38E-02  2.44E+01 -1.93E+04  2.07E+01
 
 TH10
+        6.23E+00  8.42E+01 -1.71E+01  5.31E+02 -2.91E+02 -5.67E+00  1.80E+01  9.26E+01  4.47E+01  3.31E+04
 
 TH11
+       -1.65E+01 -2.26E+01 -2.88E+01 -7.47E+01  2.98E+01  3.47E+00  3.47E+00 -1.74E+03 -4.26E-01  2.75E+01  5.92E+02
 
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
 #CPUT: Total CPU Time in Seconds,       70.453
Stop Time:
Wed Sep 29 04:52:40 CDT 2021
