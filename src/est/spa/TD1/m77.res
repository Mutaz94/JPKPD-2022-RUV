Wed Sep 29 18:28:03 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1690.14869578122        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9648E+02  3.8840E+01  1.3121E+01  7.2906E+01 -4.4924E+01  4.0465E+01  5.7581E+00 -1.2054E+00  4.5733E+01  1.3147E+01
            -1.4268E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1696.24266770164        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  1.0078E+00  1.1015E+00  1.0446E+00  9.5746E-01  1.1071E+00  1.0644E+00  1.0015E+00  1.0184E+00  7.4531E-01  9.2684E-01
             1.0848E+00
 PARAMETER:  1.0775E-01  1.9669E-01  1.4364E-01  5.6525E-02  2.0174E-01  1.6237E-01  1.0154E-01  1.1822E-01 -1.9396E-01  2.4031E-02
             1.8139E-01
 GRADIENT:  -3.1682E+01  2.4498E+01  5.7799E+00  3.0382E+01  9.8401E+00  1.2706E+01 -1.1016E+01 -9.2397E+00 -1.0063E+01 -4.6475E+00
             1.1716E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.47523956905        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0086E+00  1.0039E+00  1.1101E+00  1.0094E+00  1.0864E+00  1.0132E+00  1.2639E+00  1.2496E+00  6.5415E-01  8.3873E-01
             1.0560E+00
 PARAMETER:  1.0856E-01  1.0394E-01  2.0445E-01  1.0938E-01  1.8284E-01  1.1313E-01  3.3421E-01  3.2282E-01 -3.2442E-01 -7.5861E-02
             1.5453E-01
 GRADIENT:  -3.0485E+01  1.5665E+01 -2.3899E+00  1.0192E+01  1.7660E+01 -6.3676E+00  6.3793E+00 -2.1416E+00 -7.0171E+00 -9.1930E+00
            -6.7641E-02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.43294831696        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0231E+00  7.6000E-01  2.2641E+00  1.2058E+00  1.3279E+00  1.0297E+00  1.2873E+00  2.0949E+00  7.1307E-01  1.2114E+00
             1.0455E+00
 PARAMETER:  1.2288E-01 -1.7444E-01  9.1719E-01  2.8714E-01  3.8358E-01  1.2928E-01  3.5255E-01  8.3949E-01 -2.3817E-01  2.9176E-01
             1.4447E-01
 GRADIENT:   7.9364E+00  1.5389E+01 -3.8486E-01  3.0209E+01 -1.5379E+00  1.5396E+00  1.3137E+00 -1.4993E+00  2.1261E+00  3.6148E+00
             8.1489E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.85152675376        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0173E+00  6.4072E-01  2.8872E+00  1.2821E+00  1.3695E+00  1.0224E+00  1.3126E+00  2.4421E+00  6.9318E-01  1.2126E+00
             1.0512E+00
 PARAMETER:  1.1716E-01 -3.4516E-01  1.1603E+00  3.4849E-01  4.1445E-01  1.2218E-01  3.7203E-01  9.9287E-01 -2.6646E-01  2.9277E-01
             1.4993E-01
 GRADIENT:  -3.2741E+00  5.7651E+00  1.6882E+00  8.7368E+00 -2.9804E+00 -9.2383E-01 -9.8718E-01 -5.0455E-01 -1.3163E+00 -1.1254E+00
            -3.8854E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1700.89209419460        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0182E+00  5.2221E-01  3.8545E+00  1.3754E+00  1.4313E+00  1.0239E+00  1.3564E+00  2.8812E+00  6.8691E-01  1.2740E+00
             1.0548E+00
 PARAMETER:  1.1807E-01 -5.4969E-01  1.4492E+00  4.1877E-01  4.5860E-01  1.2359E-01  4.0486E-01  1.1582E+00 -2.7555E-01  3.4215E-01
             1.5331E-01
 GRADIENT:  -1.8546E+00  8.0493E+00  1.3897E+00  2.0990E+01 -1.6511E+00 -2.7426E-01  1.4967E-01 -2.1977E+00 -1.2445E-01  8.0649E-02
             9.6440E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1700.91368096131        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0198E+00  4.5906E-01  4.5314E+00  1.4265E+00  1.4604E+00  1.0263E+00  1.3905E+00  3.1653E+00  6.8145E-01  1.3024E+00
             1.0551E+00
 PARAMETER:  1.1965E-01 -6.7857E-01  1.6110E+00  4.5521E-01  4.7869E-01  1.2592E-01  4.2968E-01  1.2522E+00 -2.8353E-01  3.6423E-01
             1.5361E-01
 GRADIENT:   9.8645E-01  9.6036E+00  1.7305E+00  3.1520E+01 -8.8895E-01  6.1111E-01  5.1794E-01 -3.9087E+00  7.9894E-01  1.6196E+00
             2.3072E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1700.96672221610        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  1.0228E+00  3.8059E-01  5.4876E+00  1.4881E+00  1.4871E+00  1.0307E+00  1.4788E+00  3.5252E+00  6.7059E-01  1.3292E+00
             1.0547E+00
 PARAMETER:  1.2257E-01 -8.6603E-01  1.8025E+00  4.9750E-01  4.9681E-01  1.3026E-01  4.9126E-01  1.3599E+00 -2.9960E-01  3.8459E-01
             1.5329E-01
 GRADIENT:   6.7055E+00  1.0312E+01  2.2924E+00  4.2293E+01 -1.3154E+00  2.2278E+00  6.3817E-01 -5.3726E+00  1.7757E+00  3.7576E+00
             3.5558E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.28096110142        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1434
 NPARAMETR:  1.0175E+00  3.1953E-01  6.2425E+00  1.5271E+00  1.4990E+00  1.0183E+00  1.1376E+00  3.8190E+00  6.5069E-01  1.3395E+00
             1.0392E+00
 PARAMETER:  1.1732E-01 -1.0409E+00  1.9314E+00  5.2338E-01  5.0481E-01  1.1813E-01  2.2890E-01  1.4400E+00 -3.2972E-01  3.9228E-01
             1.3844E-01
 GRADIENT:  -4.4668E+00  8.8613E+00 -6.4763E-01  3.3880E+01 -7.6118E+00 -2.5938E+00 -2.2771E-01  1.9388E+00 -1.0496E+01  2.1709E+00
            -7.2149E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1701.77371198242        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1610
 NPARAMETR:  1.0176E+00  3.1940E-01  6.2423E+00  1.5264E+00  1.4991E+00  1.0207E+00  2.7284E-01  3.8190E+00  6.9581E-01  1.3394E+00
             1.0475E+00
 PARAMETER:  1.1746E-01 -1.0413E+00  1.9313E+00  5.2291E-01  5.0484E-01  1.2047E-01 -1.1989E+00  1.4400E+00 -2.6268E-01  3.9226E-01
             1.4642E-01
 GRADIENT:  -4.1962E+00  7.4864E+00 -1.4042E-01  3.6883E+01 -6.4594E+00 -1.6308E+00  5.7481E-02  5.1592E-01  2.3589E+00  2.8568E+00
            -2.0267E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1702.08648191711        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1786
 NPARAMETR:  1.0181E+00  3.1717E-01  6.2277E+00  1.5142E+00  1.4991E+00  1.0265E+00  1.0000E-02  3.8251E+00  6.8927E-01  1.3385E+00
             1.0503E+00
 PARAMETER:  1.1792E-01 -1.0483E+00  1.9290E+00  5.1488E-01  5.0488E-01  1.2615E-01 -1.4593E+01  1.4416E+00 -2.7213E-01  3.9151E-01
             1.4909E-01
 GRADIENT:  -2.7389E+00  2.4313E+00 -2.9059E-01 -1.0986E+00 -5.6069E+00  7.4752E-01  0.0000E+00  1.9644E+00 -5.5776E-01  2.6181E+00
            -1.1417E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1702.14232583047        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1952
 NPARAMETR:  1.0192E+00  3.1234E-01  6.2542E+00  1.5137E+00  1.5044E+00  1.0243E+00  1.0000E-02  3.8113E+00  6.8732E-01  1.3205E+00
             1.0538E+00
 PARAMETER:  1.1902E-01 -1.0636E+00  1.9333E+00  5.1455E-01  5.0838E-01  1.2398E-01 -1.5052E+01  1.4380E+00 -2.7495E-01  3.7801E-01
             1.5238E-01
 GRADIENT:   4.9592E+02  3.9218E+01  6.8494E+00  7.6157E+02  1.8764E+01  5.8800E+01  0.0000E+00  1.9205E+01  1.7610E+01  3.3151E+00
             1.9858E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1702.16392241184        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2128
 NPARAMETR:  1.0207E+00  2.7900E-01  6.2110E+00  1.5439E+00  1.4995E+00  1.0264E+00  1.0000E-02  3.7986E+00  6.8339E-01  1.3085E+00
             1.0514E+00
 PARAMETER:  1.2054E-01 -1.1765E+00  1.9263E+00  5.3434E-01  5.0516E-01  1.2602E-01 -1.5052E+01  1.4346E+00 -2.8069E-01  3.6886E-01
             1.5009E-01
 GRADIENT:   3.4461E+00  3.0750E+00 -6.2663E-01  8.8672E+00  4.5144E-01  6.6880E-01  0.0000E+00  9.3250E-01  2.7557E+00 -1.0217E+00
            -1.3799E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1702.16777520923        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2304
 NPARAMETR:  1.0220E+00  2.5050E-01  6.1943E+00  1.5663E+00  1.4946E+00  1.0270E+00  1.0000E-02  3.7859E+00  6.7576E-01  1.3027E+00
             1.0513E+00
 PARAMETER:  1.2176E-01 -1.2843E+00  1.9236E+00  5.4873E-01  5.0189E-01  1.2666E-01 -1.5052E+01  1.4313E+00 -2.9192E-01  3.6442E-01
             1.5003E-01
 GRADIENT:   6.5606E+00  3.7636E+00 -8.3600E-01  1.7469E+01  3.0051E-01  9.1328E-01  0.0000E+00  8.2783E-01  3.3774E+00 -1.3582E+00
            -1.6821E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1702.16867999684        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2482
 NPARAMETR:  1.0227E+00  2.2979E-01  6.1854E+00  1.5821E+00  1.4909E+00  1.0273E+00  1.0000E-02  3.7772E+00  6.6994E-01  1.2988E+00
             1.0515E+00
 PARAMETER:  1.2244E-01 -1.3706E+00  1.9222E+00  5.5877E-01  4.9941E-01  1.2695E-01 -1.5052E+01  1.4290E+00 -3.0057E-01  3.6142E-01
             1.5019E-01
 GRADIENT:   8.3762E+00  4.0142E+00 -9.2115E-01  2.2605E+01  1.6798E-01  1.0209E+00  0.0000E+00  6.7319E-01  3.5939E+00 -1.4946E+00
            -1.7598E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1702.40071559414        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2667
 NPARAMETR:  1.0193E+00  2.0884E-01  6.2120E+00  1.5822E+00  1.4886E+00  1.0251E+00  1.0000E-02  3.7736E+00  6.5920E-01  1.3059E+00
             1.0539E+00
 PARAMETER:  1.1915E-01 -1.4662E+00  1.9265E+00  5.5883E-01  4.9786E-01  1.2476E-01 -1.5052E+01  1.4280E+00 -3.1673E-01  3.6689E-01
             1.5251E-01
 GRADIENT:   2.1722E+00  1.7909E-01 -3.0617E-03 -1.8915E+01  1.4341E+00  2.8167E-01  0.0000E+00  1.6734E-02  1.5498E+00  1.8616E-01
             3.9191E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1702.46382637526        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2850
 NPARAMETR:  1.0190E+00  1.7809E-01  6.2413E+00  1.6008E+00  1.4841E+00  1.0247E+00  1.0000E-02  3.7798E+00  6.4993E-01  1.3009E+00
             1.0536E+00
 PARAMETER:  1.1884E-01 -1.6254E+00  1.9312E+00  5.7049E-01  4.9483E-01  1.2438E-01 -1.5052E+01  1.4297E+00 -3.3089E-01  3.6308E-01
             1.5224E-01
 GRADIENT:   2.1476E+00 -3.4160E-01 -1.7525E-01 -2.6997E+01  1.0349E+00  1.3649E-01  0.0000E+00  6.6793E-01  1.3677E+00 -1.3262E-02
            -1.2413E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1702.48748803540        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3035             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0207E+00  1.7238E-01  6.2365E+00  1.6081E+00  1.4811E+00  1.0256E+00  1.0000E-02  3.7815E+00  6.4593E-01  1.2988E+00
             1.0536E+00
 PARAMETER:  1.2050E-01 -1.6580E+00  1.9304E+00  5.7508E-01  4.9280E-01  1.2530E-01 -1.5052E+01  1.4301E+00 -3.3707E-01  3.6146E-01
             1.5226E-01
 GRADIENT:   5.0785E+02  1.8880E+01  5.4881E+00  9.6992E+02  1.7043E+01  5.9501E+01  0.0000E+00  1.7915E+01  2.2773E+01  2.4074E+00
             6.0496E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1702.49420724545        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     3193
 NPARAMETR:  1.0178E+00  1.6981E-01  6.2534E+00  1.6096E+00  1.4814E+00  1.0242E+00  1.0000E-02  3.7759E+00  6.4551E-01  1.2986E+00
             1.0540E+00
 PARAMETER:  1.1769E-01 -1.6731E+00  1.9331E+00  5.7601E-01  4.9297E-01  1.2387E-01 -1.5052E+01  1.4286E+00 -3.3771E-01  3.6129E-01
             1.5258E-01
 GRADIENT:  -2.3823E-01  5.0884E-01 -8.8486E-03 -1.7045E+01  2.0760E-01 -7.0652E-02  0.0000E+00  1.4851E-01  2.6532E-01  9.1506E-02
            -1.7036E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1702.50620004488        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3376             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0202E+00  1.6391E-01  6.2488E+00  1.6131E+00  1.4802E+00  1.0254E+00  1.0000E-02  3.7805E+00  6.4382E-01  1.2976E+00
             1.0541E+00
 PARAMETER:  1.1996E-01 -1.7085E+00  1.9324E+00  5.7815E-01  4.9218E-01  1.2505E-01 -1.5052E+01  1.4299E+00 -3.4034E-01  3.6051E-01
             1.5268E-01
 GRADIENT:   5.0420E+02  1.7406E+01  5.6220E+00  9.7931E+02  1.7379E+01  5.9308E+01  0.0000E+00  1.7453E+01  2.3117E+01  2.4728E+00
             9.4261E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1702.51265227056        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3554
 NPARAMETR:  1.0197E+00  1.5823E-01  6.2542E+00  1.6152E+00  1.4799E+00  1.0251E+00  1.0000E-02  3.7815E+00  6.4294E-01  1.2970E+00
             1.0542E+00
 PARAMETER:  1.1948E-01 -1.7437E+00  1.9333E+00  5.7946E-01  4.9195E-01  1.2475E-01 -1.5052E+01  1.4301E+00 -3.4170E-01  3.6005E-01
             1.5279E-01
 GRADIENT:   3.8743E+00 -4.2204E-02 -2.4664E-01 -2.4670E+01  7.2482E-02  2.8440E-01  0.0000E+00  8.7299E-01  6.6852E-01 -3.8534E-02
            -1.7977E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1702.52220763763        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     3714
 NPARAMETR:  1.0188E+00  1.5584E-01  6.2540E+00  1.6180E+00  1.4790E+00  1.0246E+00  1.0000E-02  3.7803E+00  6.4132E-01  1.2964E+00
             1.0542E+00
 PARAMETER:  1.1860E-01 -1.7589E+00  1.9332E+00  5.8118E-01  4.9134E-01  1.2426E-01 -1.5052E+01  1.4298E+00 -3.4422E-01  3.5956E-01
             1.5275E-01
 GRADIENT:   2.0109E+00  2.4336E-01 -2.4357E-01 -2.1031E+01 -2.6419E-01  9.3405E-02  0.0000E+00  8.1522E-01  1.5738E-01 -2.3001E-02
            -2.5460E-01

0ITERATION NO.:  106    OBJECTIVE VALUE:  -1702.52220763763        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     3742
 NPARAMETR:  1.0192E+00  1.5470E-01  6.2926E+00  1.6218E+00  1.4797E+00  1.0255E+00  1.0000E-02  3.7739E+00  6.4109E-01  1.2972E+00
             1.0558E+00
 PARAMETER:  1.1860E-01 -1.7589E+00  1.9332E+00  5.8118E-01  4.9134E-01  1.2426E-01 -1.5052E+01  1.4298E+00 -3.4422E-01  3.5956E-01
             1.5275E-01
 GRADIENT:  -2.8864E-01  5.3657E-02 -1.3570E-01 -3.7166E+00 -1.5676E-01 -1.1714E-01  0.0000E+00  1.7917E-01  3.9505E-02 -3.7099E-02
            -2.5375E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3742
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9924E-03 -7.1355E-05 -5.1280E-02 -1.4025E-02 -6.0121E-02
 SE:             2.9768E-02  3.0770E-05  1.7247E-02  2.8558E-02  2.0451E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4664E-01  2.0395E-02  2.9466E-03  6.2336E-01  3.2842E-03

 ETASHRINKSD(%)  2.7237E-01  9.9897E+01  4.2220E+01  4.3261E+00  3.1488E+01
 ETASHRINKVR(%)  5.4400E-01  1.0000E+02  6.6615E+01  8.4651E+00  5.3061E+01
 EBVSHRINKSD(%)  4.5766E-01  9.9900E+01  5.4982E+01  3.9562E+00  2.5270E+01
 EBVSHRINKVR(%)  9.1322E-01  1.0000E+02  7.9734E+01  7.7558E+00  4.4155E+01
 RELATIVEINF(%)  9.8178E+01  5.9743E-06  1.1837E+01  5.6760E+00  2.9637E+01
 EPSSHRINKSD(%)  4.3833E+01
 EPSSHRINKVR(%)  6.8453E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1702.5222076376317     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.37138107389353     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    59.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1702.522       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.56E-01  6.25E+00  1.62E+00  1.48E+00  1.02E+00  1.00E-02  3.78E+00  6.41E-01  1.30E+00  1.05E+00
 


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
+        1.01E+03
 
 TH 2
+       -2.15E+01  4.82E+02
 
 TH 3
+       -4.67E-01 -8.27E-01  8.95E-01
 
 TH 4
+       -1.62E+01  6.50E+02 -1.64E+00  9.67E+02
 
 TH 5
+       -6.74E-01 -4.10E+01  1.89E-01 -1.94E+01  2.22E+02
 
 TH 6
+        2.71E-02 -1.33E+00 -8.35E-02 -3.00E+00 -4.39E-01  1.87E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.18E+00 -4.19E+00 -2.87E+00 -9.28E+00 -2.19E+01 -8.31E-02  0.00E+00  1.18E+01
 
 TH 9
+        7.93E-01 -1.19E+02 -4.90E-01 -1.88E+00 -3.41E+00 -2.62E+00  0.00E+00  2.86E+00  4.23E+02
 
 TH10
+        3.76E-01 -1.64E+00  1.96E+00  2.27E+00 -1.97E+01 -8.63E-02  0.00E+00 -6.88E+00 -4.06E+00  5.53E+01
 
 TH11
+       -5.90E+00 -1.09E+01  3.72E+00 -2.35E+00  3.90E+01  2.62E+00  0.00E+00 -1.74E+01  8.61E+00  3.28E+01  2.40E+02
 
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
 #CPUT: Total CPU Time in Seconds,       67.563
Stop Time:
Wed Sep 29 18:29:12 CDT 2021
