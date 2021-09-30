Wed Sep 29 20:05:43 CDT 2021
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
$DATA ../../../../data/spa/D/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12163.3514953104        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9655E+02  2.2074E+02  8.8202E+00  2.6178E+02  2.0878E+02 -1.4673E+03 -6.2751E+02 -9.6284E+01 -7.9772E+02 -5.9502E+02
            -2.3686E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -582.201420928391        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2817E+00  1.1499E+00  9.6659E-01  1.7058E+00  1.4780E+00  2.3503E+00  1.7744E+00  1.0460E+00  1.8060E+00  1.1727E+00
             1.3416E+01
 PARAMETER:  3.4816E-01  2.3965E-01  6.6017E-02  6.3404E-01  4.9067E-01  9.5452E-01  6.7344E-01  1.4495E-01  6.9113E-01  2.5935E-01
             2.6965E+00
 GRADIENT:   2.2135E+00  1.5622E+01 -1.8658E+01  3.8545E+01 -6.2270E+00  8.6550E+01 -2.8067E+00  5.3835E+00  9.3372E+00  1.2133E+00
             1.4741E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -604.547548619295        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2430E+00  1.5616E+00  2.3944E+00  1.5186E+00  7.5096E+00  1.8957E+00  2.6690E+00  1.7068E+00  2.0808E+00  8.9943E+00
             1.2035E+01
 PARAMETER:  3.1751E-01  5.4573E-01  9.7313E-01  5.1779E-01  2.1162E+00  7.3958E-01  1.0817E+00  6.3460E-01  8.3276E-01  2.2966E+00
             2.5878E+00
 GRADIENT:   1.3261E+01  3.3430E+01  2.2711E+00  4.6779E+01  5.8071E-01  3.4872E+01  1.0468E+01  1.4698E-01  1.8086E+01 -3.5046E-01
             1.0180E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -619.502209542109        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1395E+00  1.2239E+00  2.0313E+00  1.1756E+00  6.8235E+00  1.5641E+00  1.8112E+00  3.2821E+00  1.4401E+00  9.8577E+00
             1.1598E+01
 PARAMETER:  2.3056E-01  3.0207E-01  8.0865E-01  2.6175E-01  2.0204E+00  5.4729E-01  6.9398E-01  1.2885E+00  4.6473E-01  2.3883E+00
             2.5508E+00
 GRADIENT:  -8.0800E+00 -5.3134E-01  4.1528E+00 -9.6834E+00 -2.4887E+00 -8.3330E+00  1.4347E+00  7.1601E-01  3.3058E+00  2.8795E+00
             6.8394E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -639.566016759606        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      477
 NPARAMETR:  1.0412E+00  3.9459E-01  8.0425E-01  1.5962E+00  7.8752E+00  1.6208E+00  2.8703E+00  6.5106E-01  1.2439E+00  9.0834E+00
             9.6616E+00
 PARAMETER:  1.4040E-01 -8.2990E-01 -1.1784E-01  5.6760E-01  2.1637E+00  5.8293E-01  1.1544E+00 -3.2915E-01  3.1824E-01  2.3065E+00
             2.3682E+00
 GRADIENT:  -1.4274E+01  1.6886E+01  1.6252E+00  3.9924E+01 -9.3396E+00 -3.1407E+01  5.2433E-01 -1.8307E-01 -5.3918E-01  5.7112E+00
            -5.7902E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -656.916932375536        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      652
 NPARAMETR:  9.8244E-01  1.0700E-01  4.0041E-01  1.4700E+00  7.4120E+00  1.8245E+00  1.4963E+00  1.3311E-01  9.3584E-01  7.9849E+00
             9.9421E+00
 PARAMETER:  8.2282E-02 -2.1350E+00 -8.1526E-01  4.8529E-01  2.1031E+00  7.0131E-01  5.0299E-01 -1.9166E+00  3.3693E-02  2.1776E+00
             2.3968E+00
 GRADIENT:  -1.5497E+00  7.4292E+00 -2.9203E+01  6.1404E+01 -1.2956E+01  8.3144E+00 -2.3869E-02 -1.7213E-01 -1.0348E+01  3.5637E+00
            -1.9854E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -660.130435852829        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      841             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8555E-01  6.0172E-02  4.2212E-01  1.4239E+00  7.3674E+00  1.7768E+00  1.1190E+00  1.3049E-01  8.4246E-01  8.2387E+00
             1.0240E+01
 PARAMETER:  8.5442E-02 -2.7106E+00 -7.6247E-01  4.5343E-01  2.0971E+00  6.7481E-01  2.1245E-01 -1.9364E+00 -7.1435E-02  2.2088E+00
             2.4263E+00
 GRADIENT:  -2.9147E+00  7.2853E-01  1.9607E+01 -3.9323E+01  1.2932E+01  2.5762E+01  6.0578E-03 -6.7318E-02  1.0766E+01 -1.5926E+01
             4.6598E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -660.577111376417        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      914
 NPARAMETR:  9.8422E-01  2.4670E-02  4.1855E-01  1.4237E+00  7.3850E+00  1.6689E+00  8.6044E-01  1.4351E-01  8.4089E-01  8.1315E+00
             1.0107E+01
 PARAMETER:  8.4093E-02 -3.6022E+00 -7.7097E-01  4.5327E-01  2.0994E+00  6.1216E-01 -5.0314E-02 -1.8413E+00 -7.3294E-02  2.1957E+00
             2.4133E+00
 GRADIENT:   1.5397E+01  1.6725E-01  1.5021E+01 -2.6439E+01  5.9555E+00 -7.9494E-02  1.2876E-03 -8.2823E-02  6.5772E+00  1.7014E+01
             1.6798E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -660.789769309673        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      984
 NPARAMETR:  9.7322E-01  1.0000E-02  4.1192E-01  1.4279E+00  7.2853E+00  1.6519E+00  6.2315E-01  1.9826E-01  8.3754E-01  7.9423E+00
             9.7985E+00
 PARAMETER:  7.2852E-02 -4.7275E+00 -7.8692E-01  4.5623E-01  2.0859E+00  6.0194E-01 -3.7297E-01 -1.5182E+00 -7.7290E-02  2.1722E+00
             2.3822E+00
 GRADIENT:   2.1966E+00  0.0000E+00  1.6442E+01 -2.2892E+01  1.4941E+01  5.4066E-01  9.9629E-05 -1.6887E-01  8.8004E+00 -1.7222E+01
             9.1433E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -661.177811267872        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.5616E-01  1.0000E-02  4.0141E-01  1.4362E+00  7.0222E+00  1.6719E+00  3.9788E-01  5.1051E-01  8.3182E-01  7.7440E+00
             9.4814E+00
 PARAMETER:  5.5171E-02 -5.9636E+00 -8.1278E-01  4.6202E-01  2.0491E+00  6.1394E-01 -8.2161E-01 -5.7234E-01 -8.4138E-02  2.1469E+00
             2.3493E+00
 GRADIENT:  -4.6312E-01  0.0000E+00  8.9557E+00  8.4446E+00 -3.2474E+00 -5.1702E+00  2.3001E-05 -1.9363E+00  3.2333E+00 -2.0778E+00
            -3.1729E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -661.252999771312        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  9.5577E-01  1.0000E-02  4.0113E-01  1.4364E+00  7.0159E+00  1.6719E+00  3.9423E-01  5.3650E-01  8.3168E-01  7.7412E+00
             9.4797E+00
 PARAMETER:  5.4758E-02 -5.9938E+00 -8.1347E-01  4.6215E-01  2.0482E+00  6.1398E-01 -8.3082E-01 -5.2269E-01 -8.4302E-02  2.1466E+00
             2.3492E+00
 GRADIENT:   3.8376E+01  0.0000E+00 -4.5693E+01  1.1139E+02 -1.3877E+02 -2.2377E+01 -4.0194E-05 -3.3960E+00 -1.0938E+01  8.6944E+01
            -1.7168E+02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -691.422217546242        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1294
 NPARAMETR:  9.6533E-01  1.0000E-02  3.7759E-01  1.4474E+00  6.7312E+00  1.7095E+00  4.1275E-01  2.4839E+00  6.7813E-01  7.9702E+00
             9.8210E+00
 PARAMETER:  6.4717E-02 -5.9938E+00 -8.7394E-01  4.6974E-01  2.0067E+00  6.3620E-01 -7.8491E-01  1.0098E+00 -2.8841E-01  2.1757E+00
             2.3845E+00
 GRADIENT:  -2.2327E+00  0.0000E+00  6.3054E+00  7.0951E+01  1.8216E+00 -2.5905E+01 -3.4259E-05 -1.0853E+01 -2.8693E+00  1.1737E+01
             5.1525E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -703.372923248720        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:     1400
 NPARAMETR:  9.4055E-01  1.0000E-02  3.2500E-01  1.2621E+00  6.9066E+00  1.8608E+00  4.2144E-01  2.7620E+00  7.4587E-01  7.2855E+00
             8.5201E+00
 PARAMETER:  3.8711E-02 -5.9938E+00 -1.0239E+00  3.3278E-01  2.0325E+00  7.2102E-01 -7.6408E-01  1.1160E+00 -1.9320E-01  2.0859E+00
             2.2424E+00
 GRADIENT:   1.3402E+01  0.0000E+00  7.3287E+00 -4.0711E+01  1.5859E+01 -7.1167E+00  1.5256E-04 -2.2789E+00  2.3193E+01 -3.9929E+01
            -1.2909E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -704.007436568474        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  9.4020E-01  1.0000E-02  2.9927E-01  1.2601E+00  6.9491E+00  1.8550E+00  4.2142E-01  2.7668E+00  7.4617E-01  7.2393E+00
             8.7204E+00
 PARAMETER:  3.8337E-02 -5.9938E+00 -1.1064E+00  3.3119E-01  2.0386E+00  7.1787E-01 -7.6412E-01  1.1177E+00 -1.9281E-01  2.0795E+00
             2.2657E+00
 GRADIENT:   2.6402E+01  0.0000E+00 -5.7756E+00  2.2735E+01 -4.3751E+00 -3.7242E+00 -2.7499E-05 -1.6992E+00  1.9230E+00 -9.6343E-01
            -2.1036E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -705.111780345349        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1681
 NPARAMETR:  9.0219E-01  1.0000E-02  3.0121E-01  1.2522E+00  6.9425E+00  1.8710E+00  4.2278E-01  2.7845E+00  6.8644E-01  7.2715E+00
             8.8732E+00
 PARAMETER: -2.9309E-03 -5.9938E+00 -1.0999E+00  3.2491E-01  2.0377E+00  7.2646E-01 -7.6091E-01  1.1241E+00 -2.7623E-01  2.0840E+00
             2.2830E+00
 GRADIENT:  -1.9214E+00  0.0000E+00  1.8675E+00  9.8176E+00  2.8489E+00 -4.8071E-01  4.7717E-06  4.2867E-01  5.1706E+00 -1.0184E+01
             2.5250E+00

0ITERATION NO.:   72    OBJECTIVE VALUE:  -705.191713072776        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:     1749
 NPARAMETR:  9.0171E-01  1.0000E-02  3.0275E-01  1.2501E+00  6.9979E+00  1.8644E+00  4.1956E-01  2.7995E+00  6.8357E-01  7.2272E+00
             8.9444E+00
 PARAMETER: -2.9796E-03 -5.9938E+00 -1.1000E+00  3.2481E-01  2.0373E+00  7.2633E-01 -7.6093E-01  1.1241E+00 -2.8173E-01  2.0856E+00
             2.2828E+00
 GRADIENT:   4.9760E+02  0.0000E+00 -5.4912E+01  1.6723E+02 -3.0392E+01  6.7262E+01  7.0983E-04 -4.5435E+01 -1.7111E+02  2.8478E+01
            -2.9793E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1749
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9115E-03 -1.3161E-04 -4.3296E-02 -3.5108E-02 -3.1361E-02
 SE:             2.8938E-02  3.0700E-05  1.8474E-02  1.3663E-02  7.9174E-03
 N:                     100         100         100         100         100

 P VAL.:         8.9248E-01  1.8138E-05  1.9098E-02  1.0183E-02  7.4676E-05

 ETASHRINKSD(%)  3.0550E+00  9.9897E+01  3.8110E+01  5.4228E+01  7.3476E+01
 ETASHRINKVR(%)  6.0166E+00  1.0000E+02  6.1696E+01  7.9049E+01  9.2965E+01
 EBVSHRINKSD(%)  4.3286E+00  9.9913E+01  3.1477E+01  5.5079E+01  6.8067E+01
 EBVSHRINKVR(%)  8.4699E+00  1.0000E+02  5.3046E+01  7.9821E+01  8.9803E+01
 RELATIVEINF(%)  4.7381E+01  4.8856E-06  3.4846E+00  4.7138E-01  3.7322E+00
 EPSSHRINKSD(%)  1.2600E+01
 EPSSHRINKVR(%)  2.3613E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -705.19171307277566     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       29.959113490962523     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -705.192       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.02E-01  1.00E-02  3.01E-01  1.25E+00  6.94E+00  1.87E+00  4.23E-01  2.78E+00  6.83E-01  7.28E+00  8.87E+00
 


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
+        1.58E+05
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.24E+04  0.00E+00  1.66E+04
 
 TH 4
+        4.50E+02  0.00E+00 -3.60E+03  1.10E+04
 
 TH 5
+       -1.25E+01  0.00E+00  5.53E+01 -4.49E+01  1.19E+01
 
 TH 6
+        1.39E+02  0.00E+00 -7.33E+02  5.57E+02 -1.26E+00  7.60E+02
 
 TH 7
+       -1.24E-01  0.00E+00 -9.12E-03 -7.57E-02  7.35E-04  2.55E-03  8.18E-02
 
 TH 8
+       -2.22E+02  0.00E+00  4.10E+02 -3.35E+02  4.62E+00 -7.72E+01  6.73E-03  1.82E+02
 
 TH 9
+       -8.10E+02  0.00E+00  5.11E+03 -4.15E+03  1.75E+01 -1.70E+02 -6.37E-02  5.47E+02  3.58E+04
 
 TH10
+        8.93E+00  0.00E+00 -4.84E+01  3.54E+01 -3.88E+00  1.02E+00 -9.98E-04 -4.01E+00 -1.18E+01  1.03E+01
 
 TH11
+       -1.68E+01  0.00E+00  7.95E+01 -6.74E+01  1.63E+00 -7.64E-01 -1.34E-04  7.29E+00  1.55E+01 -1.72E+00  1.03E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       40.641
Stop Time:
Wed Sep 29 20:06:25 CDT 2021
