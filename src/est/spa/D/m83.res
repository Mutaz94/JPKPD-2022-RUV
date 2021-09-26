Sat Sep 25 14:42:34 CDT 2021
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
$DATA ../../../../data/spa/D/dat83.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24803.1151104942        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.0442E+02  4.6362E+02  3.6755E+01  4.1755E+02 -7.6904E+01 -1.9632E+03 -1.0117E+03 -1.2063E+02 -1.4156E+03 -1.8216E+02
            -4.7844E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -512.277312509724        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3136E+00  9.9664E-01  8.8417E-01  1.3111E+00  1.3421E+00  1.5906E+00  1.0019E+00  9.7574E-01  8.6130E-01  9.4228E-01
             1.5116E+01
 PARAMETER:  3.7279E-01  9.6636E-02 -2.3112E-02  3.7084E-01  3.9425E-01  5.6414E-01  1.0189E-01  7.5436E-02 -4.9308E-02  4.0550E-02
             2.8158E+00
 GRADIENT:   2.6336E+01 -4.8326E+00  3.8830E-01 -2.1885E+01 -2.9615E+00  3.9029E+01  1.9385E+00  3.0131E+00  8.7116E+00  1.5398E+00
             1.0815E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -524.248782798978        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.2848E+00  6.5366E-01  8.3978E-01  1.5928E+00  1.3110E+00  1.3372E+00  8.8979E-01  2.5765E-01  3.0317E-01  5.4836E-01
             1.6190E+01
 PARAMETER:  3.5062E-01 -3.2517E-01 -7.4616E-02  5.6548E-01  3.7081E-01  3.9060E-01 -1.6769E-02 -1.2562E+00 -1.0935E+00 -5.0083E-01
             2.8844E+00
 GRADIENT:  -1.4402E+01  1.2081E+01 -3.0241E+00  2.5687E+01  3.6651E-01 -1.2116E+00  4.7453E-01  3.5498E-01  1.7684E+00  2.4306E-01
             3.5750E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -544.476340254850        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.1110E+00  1.1842E-01  2.3661E-01  1.5206E+00  2.5341E-01  1.5085E+00  5.4517E+00  1.0000E-02  1.0750E-01  1.0000E-02
             1.4106E+01
 PARAMETER:  2.0528E-01 -2.0335E+00 -1.3413E+00  5.1908E-01 -1.2728E+00  5.1113E-01  1.7959E+00 -9.9640E+00 -2.1302E+00 -5.0632E+00
             2.7466E+00
 GRADIENT:  -1.1435E+02  6.6153E+00  5.8514E+01  9.4365E+01 -7.9217E+01 -4.5560E+01  6.3670E-02  0.0000E+00  3.7719E-01  0.0000E+00
            -3.7906E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -602.427725035988        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.3357E+00  4.0713E-02  4.0842E-02  7.5365E-01  1.1412E-01  2.1105E+00  9.1110E-02  1.0000E-02  3.2573E-02  1.0000E-02
             1.3890E+01
 PARAMETER:  3.8948E-01 -3.1012E+00 -3.0980E+00 -1.8283E-01 -2.0705E+00  8.4695E-01 -2.2957E+00 -2.4296E+01 -3.3243E+00 -8.4162E+00
             2.7312E+00
 GRADIENT:   7.3664E+00 -2.6704E-01 -3.4782E-01  2.4968E+00  2.2833E+01 -3.8102E+00  4.3237E-06  0.0000E+00  3.0200E-03  0.0000E+00
             1.0507E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -604.005416725793        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.3003E+00  2.7175E-02  2.6163E-02  6.3916E-01  1.0045E-01  2.2051E+00  4.0331E-02  1.0000E-02  2.7743E-02  1.0000E-02
             1.3524E+01
 PARAMETER:  3.6259E-01 -3.5055E+00 -3.5434E+00 -3.4761E-01 -2.1981E+00  8.9078E-01 -3.1106E+00 -2.7958E+01 -3.4848E+00 -9.5700E+00
             2.7044E+00
 GRADIENT:   7.5766E-01 -9.4460E-02 -5.3386E-01  2.4780E-01  2.2438E+00 -1.0392E+00  3.9768E-07  0.0000E+00  5.9506E-04  0.0000E+00
             5.1117E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -604.006281357451        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.2990E+00  2.7378E-02  2.6092E-02  6.3715E-01  1.0027E-01  2.2129E+00  4.0086E-02  1.0000E-02  2.7888E-02  1.0000E-02
             1.3516E+01
 PARAMETER:  3.6163E-01 -3.4980E+00 -3.5461E+00 -3.5075E-01 -2.1999E+00  8.9429E-01 -3.1167E+00 -2.7967E+01 -3.4795E+00 -9.5568E+00
             2.7039E+00
 GRADIENT:   2.1410E-01 -9.4636E-02 -9.6694E-02  1.3924E-01  4.9929E-01 -1.8867E-01  4.0128E-07  0.0000E+00  6.0093E-04  0.0000E+00
             3.2616E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -604.095857274617        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  1.3045E+00  5.7649E-02  2.7132E-02  5.5834E-01  9.9984E-02  2.2493E+00  4.8063E-02  1.0000E-02  3.3187E-02  1.0000E-02
             1.3471E+01
 PARAMETER:  3.6582E-01 -2.7534E+00 -3.5071E+00 -4.8278E-01 -2.2027E+00  9.1061E-01 -2.9352E+00 -2.6460E+01 -3.3056E+00 -7.7253E+00
             2.7005E+00
 GRADIENT:  -4.4442E+00 -8.2636E-01  5.3240E+00 -1.8673E+00 -1.4535E+01  6.9571E+00  1.2790E-04  0.0000E+00  1.3887E-03  0.0000E+00
            -4.7845E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -604.838255441952        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      596
 NPARAMETR:  1.3255E+00  8.9094E-02  2.2883E-02  5.3869E-01  9.8926E-02  2.1905E+00  1.0728E-01  1.0000E-02  2.2230E-02  1.0000E-02
             1.3468E+01
 PARAMETER:  3.8180E-01 -2.3181E+00 -3.6774E+00 -5.1862E-01 -2.2134E+00  8.8413E-01 -2.1323E+00 -2.5750E+01 -3.7063E+00 -6.0161E+00
             2.7003E+00
 GRADIENT:   3.4541E+00 -5.9935E-01 -4.1206E-02 -2.3568E+00  5.0874E+00 -1.9400E+00  6.0242E-03  0.0000E+00  3.3996E-04  0.0000E+00
            -2.1280E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -605.224937081636        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.3068E+00  9.4810E-02  2.2235E-02  6.7478E-01  9.8035E-02  2.2508E+00  4.9180E-01  1.0000E-02  1.3519E-02  1.0000E-02
             1.3556E+01
 PARAMETER:  3.6756E-01 -2.2559E+00 -3.7061E+00 -2.9336E-01 -2.2224E+00  9.1128E-01 -6.0968E-01 -2.4012E+01 -4.2036E+00 -4.8781E+00
             2.7069E+00
 GRADIENT:   4.8607E-01 -2.4062E-01  5.5724E-01 -1.3148E+00 -9.8966E-01  4.3693E-01  1.7482E-01  0.0000E+00  4.6562E-05  0.0000E+00
             2.9051E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -605.785001069188        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      741
 NPARAMETR:  1.2505E+00  1.0081E-01  2.3982E-02  1.1810E+00  9.7302E-02  2.3598E+00  1.5101E-01  1.0000E-02  9.7046E-02  4.7039E-02
             1.3491E+01
 PARAMETER:  3.2357E-01 -2.1945E+00 -3.6305E+00  2.6635E-01 -2.2299E+00  9.5858E-01 -1.7904E+00 -2.0573E+01 -2.2326E+00 -2.9568E+00
             2.7020E+00
 GRADIENT:  -4.5217E-01 -6.4693E-01  1.1101E+00 -1.0465E+00 -1.1295E+00 -7.8572E-02  2.2664E-02  0.0000E+00  6.8159E-06  1.6091E-01
            -6.4140E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -606.289230515036        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      814
 NPARAMETR:  1.1983E+00  1.0460E-01  2.7236E-02  2.1926E+00  9.5601E-02  2.5050E+00  2.2026E-02  1.0000E-02  1.5032E+00  3.4997E-02
             1.3547E+01
 PARAMETER:  2.8092E-01 -2.1577E+00 -3.5032E+00  8.8508E-01 -2.2476E+00  1.0183E+00 -3.7155E+00 -1.7200E+01  5.0758E-01 -3.2525E+00
             2.7062E+00
 GRADIENT:   2.0263E+00 -3.1339E-01 -7.0787E-01 -2.6178E-01  1.9749E+00  9.1675E-01  5.6029E-04  0.0000E+00  1.4633E-01  9.0029E-02
             1.3302E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -606.649415288657        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.1076E+00  1.0279E-01  3.3376E-02  4.1474E+00  9.4090E-02  2.6960E+00  3.7753E-01  1.0000E-02  1.9744E+00  1.0000E-02
             1.3563E+01
 PARAMETER:  2.0222E-01 -2.1751E+00 -3.2999E+00  1.5225E+00 -2.2635E+00  1.0918E+00 -8.7411E-01 -1.2804E+01  7.8028E-01 -5.3093E+00
             2.7073E+00
 GRADIENT:   5.2410E-01 -8.8055E-03 -3.0265E-01 -4.4725E-01  4.7043E-01  8.5673E-01  1.5339E-01  0.0000E+00  2.8610E-01  0.0000E+00
             1.7332E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -606.755629635067        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      957
 NPARAMETR:  1.0813E+00  1.0288E-01  3.5521E-02  4.9119E+00  9.3785E-02  2.7573E+00  3.1638E-01  1.0000E-02  2.0779E+00  2.6872E-02
             1.3577E+01
 PARAMETER:  1.7815E-01 -2.1742E+00 -3.2376E+00  1.6917E+00 -2.2667E+00  1.1142E+00 -1.0508E+00 -1.1349E+01  8.3137E-01 -3.5167E+00
             2.7084E+00
 GRADIENT:   6.3842E-01  1.3851E-01  3.3329E-01 -5.2910E-01 -6.5853E-01  4.4043E-01  1.0534E-01  0.0000E+00  2.5046E-01  5.4784E-02
             1.6779E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -606.938892980857        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:     1034
 NPARAMETR:  9.9577E-01  9.2895E-02  4.1322E-02  8.4269E+00  9.1461E-02  2.9316E+00  1.2203E-01  1.0000E-02  2.8446E+00  2.4591E-02
             1.3590E+01
 PARAMETER:  9.5764E-02 -2.2763E+00 -3.0864E+00  2.2314E+00 -2.2918E+00  1.1756E+00 -2.0035E+00 -7.8414E+00  1.1454E+00 -3.6054E+00
             2.7093E+00
 GRADIENT:  -2.0274E-01 -6.8249E-01  9.8716E-01 -7.1794E-01 -8.9350E-01 -6.3388E-01  1.0541E-02  0.0000E+00  9.3226E-01  4.6802E-02
             2.1795E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -607.114761390709        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1149
 NPARAMETR:  9.5140E-01  9.7668E-02  4.5401E-02  1.0687E+01  9.1608E-02  3.0903E+00  3.8981E-02  1.0000E-02  2.8303E+00  1.0000E-02
             1.3609E+01
 PARAMETER:  5.0176E-02 -2.2262E+00 -2.9922E+00  2.4690E+00 -2.2902E+00  1.2283E+00 -3.1447E+00 -6.1007E+00  1.1404E+00 -4.9784E+00
             2.7107E+00
 GRADIENT:   1.8633E+00  1.2351E-01  2.7221E-01  3.9998E-01 -3.7171E-01 -1.0978E+00  1.2045E-03  0.0000E+00 -3.1225E-03  0.0000E+00
            -2.4992E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -607.139694260103        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1324
 NPARAMETR:  9.2526E-01  9.7189E-02  4.6972E-02  1.1337E+01  9.1717E-02  3.1597E+00  1.9476E-02  1.0000E-02  2.8094E+00  1.0000E-02
             1.3642E+01
 PARAMETER:  2.2317E-02 -2.2311E+00 -2.9582E+00  2.5280E+00 -2.2890E+00  1.2505E+00 -3.8386E+00 -5.6611E+00  1.1330E+00 -5.4245E+00
             2.7131E+00
 GRADIENT:  -4.0671E-03  1.6569E-03  2.4660E-03 -2.3029E-03 -1.0187E-03  5.5954E-03  2.9098E-04  0.0000E+00  1.8599E-03  0.0000E+00
             2.7784E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -607.139801534582        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1503
 NPARAMETR:  9.2541E-01  9.7187E-02  4.6960E-02  1.1330E+01  9.1717E-02  3.1592E+00  1.0000E-02  1.0000E-02  2.8089E+00  1.0000E-02
             1.3642E+01
 PARAMETER:  2.2485E-02 -2.2311E+00 -2.9585E+00  2.5275E+00 -2.2890E+00  1.2503E+00 -4.6568E+00 -5.7206E+00  1.1328E+00 -5.6144E+00
             2.7131E+00
 GRADIENT:  -8.4602E-04 -1.0130E-04 -1.1244E-03  3.0187E-04  4.3602E-04 -3.3588E-04  0.0000E+00  0.0000E+00 -2.1969E-04  0.0000E+00
             1.3466E-03

0ITERATION NO.:   86    OBJECTIVE VALUE:  -607.139801534582        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1525
 NPARAMETR:  9.2541E-01  9.7187E-02  4.6960E-02  1.1330E+01  9.1717E-02  3.1592E+00  1.0000E-02  1.0000E-02  2.8089E+00  1.0000E-02
             1.3642E+01
 PARAMETER:  2.2485E-02 -2.2311E+00 -2.9585E+00  2.5275E+00 -2.2890E+00  1.2503E+00 -4.6568E+00 -5.7206E+00  1.1328E+00 -5.6144E+00
             2.7131E+00
 GRADIENT:  -8.4602E-04 -1.0130E-04 -1.1244E-03  3.0187E-04  4.3602E-04 -3.3588E-04  0.0000E+00  0.0000E+00 -2.1969E-04  0.0000E+00
             1.3466E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1525
 NO. OF SIG. DIGITS IN FINAL EST.:  4.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2716E-02 -2.3258E-05 -1.4052E-04 -2.8857E-02  2.9605E-04
 SE:             2.6884E-02  1.5067E-05  6.7563E-05  1.1149E-02  2.3679E-04
 N:                     100         100         100         100         100

 P VAL.:         6.3623E-01  1.2268E-01  3.7542E-02  9.6454E-03  2.1121E-01

 ETASHRINKSD(%)  9.9354E+00  9.9950E+01  9.9774E+01  6.2649E+01  9.9207E+01
 ETASHRINKVR(%)  1.8884E+01  1.0000E+02  9.9999E+01  8.6049E+01  9.9994E+01
 EBVSHRINKSD(%)  7.8524E+00  9.9921E+01  9.9735E+01  6.6105E+01  9.9002E+01
 EBVSHRINKVR(%)  1.5088E+01  1.0000E+02  9.9999E+01  8.8511E+01  9.9990E+01
 RELATIVEINF(%)  3.2363E+01  5.1442E-05  1.6221E-04  2.6740E+00  3.3059E-03
 EPSSHRINKSD(%)  4.6708E+00
 EPSSHRINKVR(%)  9.1234E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -607.13980153458238     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       128.01102502915580     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.58
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -607.140       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.25E-01  9.72E-02  4.70E-02  1.13E+01  9.17E-02  3.16E+00  1.00E-02  1.00E-02  2.81E+00  1.00E-02  1.36E+01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.13E-01
 
 TH 2
+        6.08E+00  3.28E+02
 
 TH 3
+        5.51E+01  2.97E+03  2.70E+04
 
 TH 4
+       -8.85E-02 -4.77E+00 -4.33E+01  6.94E-02
 
 TH 5
+       -1.04E+02 -5.63E+03 -5.11E+04  8.19E+01  9.67E+04
 
 TH 6
+       -1.60E-01 -8.62E+00 -7.82E+01  1.25E-01  1.48E+02  2.27E-01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.40E-02  1.83E+00  1.66E+01 -2.67E-02 -3.15E+01 -4.82E-02  0.00E+00  0.00E+00  1.03E-02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.59E-02 -4.63E+00 -4.20E+01  6.74E-02  7.95E+01  1.22E-01  0.00E+00  0.00E+00 -2.59E-02  0.00E+00  6.55E-02
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.09E+02
 
 TH 2
+       -3.15E+00  1.16E+03
 
 TH 3
+        6.32E+02  2.02E+03  4.05E+04
 
 TH 4
+        1.27E+00 -2.66E+00 -4.15E+01  1.41E-01
 
 TH 5
+        1.53E+02 -5.17E+03 -3.79E+04  7.06E+01  8.80E+04
 
 TH 6
+        2.17E+00 -9.44E+00 -2.05E+02 -3.53E-01  7.29E+01  1.11E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.12E+00  2.60E+00  3.95E+01 -2.72E-01 -1.72E+01  1.22E+00  0.00E+00  0.00E+00  1.31E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.66E+00 -7.27E+00 -8.32E+01 -3.56E-02  5.04E+01  3.14E-01  0.00E+00  0.00E+00  1.55E-01  0.00E+00  2.12E+00
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.11E+02
 
 TH 2
+       -6.08E+00  3.04E+02
 
 TH 3
+        4.70E+02  1.53E+03  2.08E+04
 
 TH 4
+        1.56E+00 -2.67E+00 -1.72E+01  1.38E-01
 
 TH 5
+        3.61E+02 -2.82E+03 -3.12E+04  7.07E+01  8.19E+04
 
 TH 6
+        8.16E+00 -4.88E+00 -1.28E+02 -5.76E-01  4.45E+01  1.56E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.44E+00 -3.84E-01 -3.31E+01 -2.26E-01  1.56E+01  1.58E+00  0.00E+00  0.00E+00  1.07E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.42E+00 -1.70E+00 -8.93E+01 -1.28E-01  4.53E+01  3.74E+00  0.00E+00  0.00E+00  1.68E-01  0.00E+00  1.30E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.829
Stop Time:
Sat Sep 25 14:43:05 CDT 2021
