Sat Sep 25 01:07:34 CDT 2021
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
$DATA ../../../../data/int/SL2/dat30.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      998
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

 TOT. NO. OF OBS RECS:      898
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2023.57698992349        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6924E+01  4.0825E+01  1.6124E+02  1.0517E+02  1.2021E+02  4.1691E+00 -7.0194E+01 -2.5362E+02 -9.0988E+01 -2.2098E+01
            -3.3212E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3032.45457814719        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0079E+00  1.3489E+00  9.6001E-01  7.5877E-01  1.2379E+00  9.9239E-01  1.0454E+00  1.0667E+00  1.0598E+00  8.9320E-01
             1.9838E+00
 PARAMETER:  1.0782E-01  3.9930E-01  5.9193E-02 -1.7606E-01  3.1340E-01  9.2359E-02  1.4438E-01  1.6453E-01  1.5811E-01 -1.2943E-02
             7.8503E-01
 GRADIENT:  -2.1166E+01 -2.0884E+01 -3.3580E+01 -1.5350E+00  4.6035E+01 -2.5752E+00  9.4970E+00 -3.5165E+00 -1.4550E+01 -3.4628E+01
            -1.8256E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3043.43044066171        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0062E+00  1.7744E+00  1.1353E+00  5.7688E-01  1.5648E+00  9.9244E-01  8.0684E-01  1.0094E+00  1.4399E+00  1.1669E+00
             1.9900E+00
 PARAMETER:  1.0617E-01  6.7344E-01  2.2692E-01 -4.5012E-01  5.4776E-01  9.2409E-02 -1.1463E-01  1.0931E-01  4.6459E-01  2.5435E-01
             7.8813E-01
 GRADIENT:  -2.6718E+01  1.2822E+02 -5.0615E+00  6.0469E+01  1.9534E+01 -3.1190E+00 -7.3727E+00 -7.2108E+00  1.1683E+00 -2.7506E+01
            -1.6467E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3054.90307673617        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0147E+00  1.5515E+00  1.4537E+00  6.7322E-01  1.5049E+00  9.9514E-01  8.3206E-01  2.7997E+00  1.2641E+00  1.1921E+00
             2.0258E+00
 PARAMETER:  1.1456E-01  5.3919E-01  4.7410E-01 -2.9568E-01  5.0870E-01  9.5133E-02 -8.3847E-02  1.1295E+00  3.3438E-01  2.7573E-01
             8.0598E-01
 GRADIENT:  -7.3062E+00  8.1375E+00 -2.4994E+01  2.2992E+01  1.9240E+01 -1.0301E+00  8.4524E-01  2.3172E+00  7.6118E+00 -3.9998E+00
            -6.6801E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3057.36397483624        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0158E+00  1.5259E+00  1.5733E+00  6.8762E-01  1.5048E+00  9.9572E-01  8.4788E-01  3.0464E+00  1.1976E+00  1.1856E+00
             2.0425E+00
 PARAMETER:  1.1568E-01  5.2260E-01  5.5318E-01 -2.7452E-01  5.0869E-01  9.5712E-02 -6.5018E-02  1.2140E+00  2.8031E-01  2.7021E-01
             8.1417E-01
 GRADIENT:  -4.9786E+00  2.6615E+00 -3.1834E+01  1.6469E+01  8.4627E+00 -7.6416E-01  1.5274E-01  1.9810E+00  5.2517E+00  6.8469E-01
            -4.5161E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3063.95334684992        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  1.0170E+00  1.5260E+00  2.4171E+00  7.0341E-01  1.5053E+00  9.9707E-01  9.0208E-01  3.0444E+00  1.1029E+00  1.2034E+00
             2.0434E+00
 PARAMETER:  1.1687E-01  5.2262E-01  9.8256E-01 -2.5182E-01  5.0896E-01  9.7062E-02 -3.0540E-03  1.2133E+00  1.9798E-01  2.8517E-01
             8.1462E-01
 GRADIENT:  -1.6744E+00  5.3215E+01 -1.3917E+00 -1.4955E+00 -5.2412E+01 -2.9338E-02  2.7395E+00 -9.0722E+00  1.3652E+00 -6.1960E-01
            -4.1195E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3064.09457174013        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0248E+00  1.5260E+00  2.4875E+00  7.0961E-01  1.5053E+00  1.0016E+00  8.8981E-01  3.0444E+00  1.0957E+00  1.2114E+00
             2.0434E+00
 PARAMETER:  1.2453E-01  5.2262E-01  1.0113E+00 -2.4304E-01  5.0897E-01  1.0160E-01 -1.6744E-02  1.2133E+00  1.9140E-01  2.9176E-01
             8.1462E-01
 GRADIENT:   1.3254E-02  4.3876E+01 -5.7773E-03  3.3323E-02 -5.9211E+01  9.9205E-02  1.7925E-03 -1.0647E+01  1.9145E-02  6.4605E-02
            -4.3007E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3065.62560128379        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      746
 NPARAMETR:  1.0233E+00  1.4806E+00  2.4867E+00  7.0966E-01  1.5051E+00  9.9912E-01  8.8493E-01  3.5102E+00  1.0743E+00  1.2105E+00
             2.0775E+00
 PARAMETER:  1.2299E-01  4.9241E-01  1.0110E+00 -2.4297E-01  5.0884E-01  9.9122E-02 -2.2244E-02  1.3557E+00  1.7168E-01  2.9100E-01
             8.3115E-01
 GRADIENT:   1.1695E+01  6.0351E+00 -5.6194E+00 -2.5137E+01 -4.5922E+01  1.0519E+00 -4.0164E-01  1.6254E+00  1.0112E+00  3.0089E+00
             5.6181E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3065.68676013059        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      924
 NPARAMETR:  1.0250E+00  1.4908E+00  2.4867E+00  7.0966E-01  1.5051E+00  1.0008E+00  9.0133E-01  3.4898E+00  1.0477E+00  1.2105E+00
             2.0752E+00
 PARAMETER:  1.2474E-01  4.9930E-01  1.0110E+00 -2.4297E-01  5.0884E-01  1.0082E-01 -3.8861E-03  1.3499E+00  1.4656E-01  2.9100E-01
             8.3003E-01
 GRADIENT:   1.7668E-01  1.4849E-01 -5.4968E+00 -2.3269E+01 -5.3177E+01  8.3257E-03  6.7441E-02 -4.8877E-02  4.4019E-02  2.1531E+00
            -7.7943E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3068.57720468707        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  1.0208E+00  1.4733E+00  2.8032E+00  7.3154E-01  1.5784E+00  9.9854E-01  9.0364E-01  3.5462E+00  1.0316E+00  1.2007E+00
             2.0741E+00
 PARAMETER:  1.2055E-01  4.8748E-01  1.1308E+00 -2.1260E-01  5.5642E-01  9.8540E-02 -1.3279E-03  1.3659E+00  1.3110E-01  2.8288E-01
             8.2952E-01
 GRADIENT:   5.2681E+00  1.9306E+01 -5.7048E+00 -6.6076E+00 -4.9296E+00  5.2896E-01 -7.2045E-01 -2.1198E+00  9.4141E-01 -8.2659E-01
             2.7355E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3068.61142228379        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  1.0207E+00  1.4723E+00  2.8104E+00  7.3232E-01  1.5789E+00  9.9852E-01  9.0404E-01  3.5494E+00  1.0304E+00  1.2005E+00
             2.0740E+00
 PARAMETER:  1.2049E-01  4.8684E-01  1.1333E+00 -2.1154E-01  5.5673E-01  9.8514E-02 -8.8161E-04  1.3668E+00  1.2991E-01  2.8273E-01
             8.2948E-01
 GRADIENT:   5.1166E+00  1.9341E+01 -5.6111E+00 -6.4885E+00 -4.4588E+00  5.1586E-01 -7.5459E-01 -2.1063E+00  9.1441E-01 -8.1270E-01
             2.7233E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3068.61490617323        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  1.0207E+00  1.4721E+00  2.8104E+00  7.3244E-01  1.5789E+00  1.0018E+00  9.0436E-01  3.5482E+00  1.0301E+00  1.2006E+00
             2.0744E+00
 PARAMETER:  1.2052E-01  4.8672E-01  1.1333E+00 -2.1137E-01  5.5672E-01  1.0185E-01 -5.2600E-04  1.3664E+00  1.2963E-01  2.8282E-01
             8.2969E-01
 GRADIENT:  -9.6840E+00  3.0474E+00 -6.1623E+00 -9.9052E+00 -8.9684E+00  1.2617E-01 -9.6989E-01 -2.7770E+00  7.4071E-01 -1.0964E+00
             1.7065E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3068.62866514295        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1488
 NPARAMETR:  1.0207E+00  1.4718E+00  2.8147E+00  7.3266E-01  1.5789E+00  1.0014E+00  9.0460E-01  3.5482E+00  1.0297E+00  1.2007E+00
             2.0744E+00
 PARAMETER:  1.2052E-01  4.8649E-01  1.1349E+00 -2.1107E-01  5.5673E-01  1.0139E-01 -2.5848E-04  1.3664E+00  1.2928E-01  2.8291E-01
             8.2966E-01
 GRADIENT:  -9.6786E+00  3.0754E+00 -6.0528E+00 -9.9780E+00 -8.9946E+00 -4.6299E-02 -9.6301E-01 -2.8158E+00  7.3383E-01 -1.0606E+00
             1.6602E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3068.64495029816        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1658
 NPARAMETR:  1.0226E+00  1.4717E+00  2.8146E+00  7.3265E-01  1.5789E+00  1.0012E+00  9.0738E-01  3.5476E+00  1.0297E+00  1.2019E+00
             2.0746E+00
 PARAMETER:  1.2231E-01  4.8643E-01  1.1348E+00 -2.1108E-01  5.5674E-01  1.0121E-01  2.8089E-03  1.3663E+00  1.2929E-01  2.8387E-01
             8.2977E-01
 GRADIENT:  -5.7996E+00  2.8896E+00 -6.0541E+00 -1.0060E+01 -9.0510E+00 -7.3910E-02 -5.3796E-01 -2.8068E+00  9.6733E-01 -8.4783E-01
             1.8740E+00

0ITERATION NO.:   69    OBJECTIVE VALUE:  -3068.65023031564        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1793
 NPARAMETR:  1.0235E+00  1.4719E+00  2.8156E+00  7.3261E-01  1.5787E+00  1.0013E+00  9.0868E-01  3.5475E+00  1.0297E+00  1.2018E+00
             2.0746E+00
 PARAMETER:  1.2328E-01  4.8643E-01  1.1348E+00 -2.1108E-01  5.5674E-01  1.0125E-01  4.2708E-03  1.3663E+00  1.2929E-01  2.8387E-01
             8.2977E-01
 GRADIENT:   1.6783E+05 -4.2542E+04 -9.1202E+03  9.8011E+04  3.7154E+04 -3.4250E-02  2.0691E+05  7.4168E+02  8.0018E+04  7.2888E+04
            -1.3421E+03
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.8         3.3         4.9         3.3         3.3
                    4.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1793
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7371E-03 -2.5577E-02 -2.2559E-02  2.4506E-02 -2.9107E-02
 SE:             2.9645E-02  2.2579E-02  1.8202E-02  2.1178E-02  2.3190E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2644E-01  2.5732E-01  2.1520E-01  2.4721E-01  2.0942E-01

 ETASHRINKSD(%)  6.8678E-01  2.4356E+01  3.9022E+01  2.9051E+01  2.2311E+01
 ETASHRINKVR(%)  1.3688E+00  4.2780E+01  6.2817E+01  4.9662E+01  3.9644E+01
 EBVSHRINKSD(%)  9.8944E-01  2.4593E+01  4.3023E+01  3.2997E+01  1.8374E+01
 EBVSHRINKVR(%)  1.9691E+00  4.3137E+01  6.7537E+01  5.5106E+01  3.3372E+01
 RELATIVEINF(%)  9.8001E+01  7.7451E+00  1.4889E+01  6.0357E+00  3.6083E+01
 EPSSHRINKSD(%)  1.8143E+01
 EPSSHRINKVR(%)  3.2994E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3068.6502303156367     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1418.2366246800445     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    53.83
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3068.650       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.47E+00  2.81E+00  7.33E-01  1.58E+00  1.00E+00  9.09E-01  3.55E+00  1.03E+00  1.20E+00  2.07E+00
 


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
+        3.25E+08
 
 TH 2
+       -3.20E+02  2.07E+07
 
 TH 3
+       -6.96E+01  8.21E+02  5.07E+05
 
 TH 4
+       -1.90E+04  4.92E+07 -1.10E+07  4.44E+08
 
 TH 5
+        2.50E+02  8.65E+06  9.61E+02  4.01E+07  1.37E+07
 
 TH 6
+        1.22E+03 -2.16E+02 -4.78E+01  9.77E+02  1.73E+02  1.89E+02
 
 TH 7
+        2.25E+03 -3.88E+02 -8.83E+01  1.83E+03  3.22E+02  1.68E+03  6.26E+08
 
 TH 8
+       -4.59E+00 -4.94E+00 -6.50E+00 -1.66E+01  1.21E+06  1.71E+01  4.09E+01  6.71E+05
 
 TH 9
+       -3.53E+03  6.25E+02  1.39E+02  2.51E+08 -5.01E+02  1.15E+03  2.17E+03 -1.21E+01  2.92E+08
 
 TH10
+       -7.63E+03  1.33E+03  3.01E+02  1.03E+08 -1.11E+03  4.47E+02  1.67E+08 -1.48E+02 -7.23E+03  4.44E+07
 
 TH11
+        1.01E+00  3.57E+03 -1.99E+01  5.67E+01 -1.16E+03 -4.72E+01 -1.11E+02 -6.19E+05  4.89E+01  4.25E+02  3.58E+06
 
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
 #CPUT: Total CPU Time in Seconds,       69.896
Stop Time:
Sat Sep 25 01:08:45 CDT 2021
