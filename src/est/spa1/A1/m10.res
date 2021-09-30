Wed Sep 29 21:49:50 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat10.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1840.75613364818        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7150E+02  3.7307E+01 -1.0579E+01  7.1544E+01  5.7996E+01  2.3742E+01  3.5245E+00  7.3328E+00 -6.4427E+00  7.8072E+00
            -4.7988E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1905.42289544483        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.3805E-01  9.1829E-01  1.0525E+00  1.0666E+00  9.6660E-01  1.0546E+00  8.9052E-01  8.6161E-01  1.0556E+00  7.2187E-01
             1.8965E+00
 PARAMETER:  3.6052E-02  1.4758E-02  1.5116E-01  1.6446E-01  6.6028E-02  1.5314E-01 -1.5949E-02 -4.8958E-02  1.5414E-01 -2.2591E-01
             7.4001E-01
 GRADIENT:   1.5215E+01  2.5299E+00 -1.4644E+01  4.0948E+01  3.3352E+01  1.9929E+01  5.5540E+00  4.1440E+00  1.1447E+01  1.1721E+00
             1.6741E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1909.66909079040        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  9.4009E-01  1.0171E+00  6.4896E-01  9.9584E-01  7.5537E-01  1.0329E+00  8.4687E-01  4.5781E-01  1.0648E+00  5.4754E-01
             1.8185E+00
 PARAMETER:  3.8223E-02  1.1693E-01 -3.3238E-01  9.5828E-02 -1.8055E-01  1.3235E-01 -6.6206E-02 -6.8130E-01  1.6275E-01 -5.0231E-01
             6.9802E-01
 GRADIENT:   1.5826E+01  5.8614E+01  1.4006E+01  5.3094E+01 -3.9843E+01  6.5940E+00 -8.0147E+00  1.0242E+00  9.3611E+00 -1.8000E+00
             1.4298E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1921.29154993035        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      248
 NPARAMETR:  9.1458E-01  8.6901E-01  5.0177E-01  1.0121E+00  6.2329E-01  1.0222E+00  1.1350E+00  1.5990E-01  9.4924E-01  5.6843E-01
             1.5238E+00
 PARAMETER:  1.0711E-02 -4.0400E-02 -5.8961E-01  1.1205E-01 -3.7275E-01  1.2200E-01  2.2663E-01 -1.7332E+00  4.7908E-02 -4.6487E-01
             5.2123E-01
 GRADIENT:  -1.6561E+02 -9.0273E+00 -3.3926E+01 -2.1161E+01  1.3503E+01 -3.3376E+01  1.3661E+00  4.8355E-01  2.8076E+00  7.4216E+00
             3.4657E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1935.34680789165        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  9.7933E-01  5.8749E-01  7.6850E-01  1.2314E+00  6.8309E-01  1.0725E+00  1.4800E+00  9.6946E-02  8.7329E-01  7.3516E-01
             1.4564E+00
 PARAMETER:  7.9115E-02 -4.3189E-01 -1.6332E-01  3.0814E-01 -2.8113E-01  1.7002E-01  4.9201E-01 -2.2336E+00 -3.5483E-02 -2.0767E-01
             4.7595E-01
 GRADIENT:  -7.4108E+00  6.4628E+00 -3.0231E-01  1.1617E+01  3.9308E+00  1.4544E+00  3.3251E-01 -3.7100E-02  2.8791E+00 -4.6379E-01
            -6.8641E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1936.47928764726        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.8033E-01  3.4849E-01  7.6854E-01  1.3549E+00  6.1603E-01  1.0752E+00  2.0756E+00  2.7498E-02  8.0586E-01  7.4855E-01
             1.4688E+00
 PARAMETER:  8.0129E-02 -9.5414E-01 -1.6326E-01  4.0376E-01 -3.8446E-01  1.7248E-01  8.3024E-01 -3.4937E+00 -1.1584E-01 -1.8962E-01
             4.8444E-01
 GRADIENT:   3.7892E+00  2.0461E+00 -4.7034E-01  5.7362E+00 -1.9169E+00  4.3051E+00  2.7577E-01 -6.2881E-03  2.4102E-01 -7.4115E-01
            -1.0067E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1936.56386290925        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      776
 NPARAMETR:  9.7572E-01  2.6701E-01  8.0376E-01  1.4036E+00  6.1655E-01  1.0630E+00  2.4257E+00  1.4030E-02  7.8760E-01  7.7734E-01
             1.4740E+00
 PARAMETER:  7.5422E-02 -1.2205E+00 -1.1846E-01  4.3902E-01 -3.8361E-01  1.6109E-01  9.8612E-01 -4.1666E+00 -1.3876E-01 -1.5188E-01
             4.8799E-01
 GRADIENT:  -7.3191E-01  1.2625E+00  2.0837E+00  4.3686E+00 -3.3071E+00  6.9807E-01  8.4524E-02 -1.9714E-03 -7.0020E-02  2.4303E-01
             4.9227E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1936.57962612571        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      964             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7628E-01  2.6460E-01  8.0273E-01  1.4010E+00  6.1671E-01  1.0616E+00  2.4390E+00  1.5878E-02  7.8728E-01  7.7638E-01
             1.4732E+00
 PARAMETER:  7.5995E-02 -1.2295E+00 -1.1974E-01  4.3716E-01 -3.8336E-01  1.5977E-01  9.9159E-01 -4.0428E+00 -1.3917E-01 -1.5312E-01
             4.8740E-01
 GRADIENT:   1.8784E+02  2.0218E+01  2.8230E+00  3.2725E+02  2.4700E+01  3.8450E+01  8.0182E+00 -2.9011E-04  6.1357E+00  6.6049E-01
             4.9833E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1936.58065058067        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  9.7632E-01  2.6493E-01  8.0233E-01  1.4006E+00  6.1649E-01  1.0615E+00  2.4400E+00  1.7911E-02  7.8698E-01  7.7607E-01
             1.4730E+00
 PARAMETER:  7.6038E-02 -1.2283E+00 -1.2024E-01  4.3689E-01 -3.8371E-01  1.5966E-01  9.9201E-01 -3.9223E+00 -1.3956E-01 -1.5352E-01
             4.8728E-01
 GRADIENT:   6.4311E-01  6.3256E-02  7.4434E-02 -3.5945E+00  7.2073E-01  1.7406E-01  1.0386E-01 -3.1885E-03  1.6726E-02  4.0078E-02
             1.1803E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1936.58201715367        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1339             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7633E-01  2.6584E-01  8.0190E-01  1.4002E+00  6.1600E-01  1.0615E+00  2.4395E+00  2.2558E-02  7.8705E-01  7.7557E-01
             1.4729E+00
 PARAMETER:  7.6047E-02 -1.2249E+00 -1.2077E-01  4.3662E-01 -3.8451E-01  1.5966E-01  9.9180E-01 -3.6917E+00 -1.3946E-01 -1.5415E-01
             4.8723E-01
 GRADIENT:   1.8792E+02  2.0594E+01  4.2407E+00  3.2671E+02  2.2741E+01  3.8382E+01  8.2422E+00 -1.0235E-03  6.0411E+00  6.7013E-01
             4.8353E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1936.58331631057        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1521             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7634E-01  2.6520E-01  8.0085E-01  1.4001E+00  6.1603E-01  1.0615E+00  2.4365E+00  2.7051E-02  7.8706E-01  7.7524E-01
             1.4729E+00
 PARAMETER:  7.6058E-02 -1.2273E+00 -1.2209E-01  4.3654E-01 -3.8446E-01  1.5967E-01  9.9057E-01 -3.5100E+00 -1.3944E-01 -1.5458E-01
             4.8724E-01
 GRADIENT:   1.8794E+02  2.0157E+01  2.1912E+00  3.2627E+02  2.5597E+01  3.8377E+01  8.0475E+00 -1.6250E-03  6.0248E+00  6.5366E-01
             4.9139E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1936.58544137161        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1706             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7633E-01  2.6568E-01  8.0107E-01  1.4002E+00  6.1560E-01  1.0615E+00  2.4396E+00  3.2150E-02  7.8711E-01  7.7505E-01
             1.4728E+00
 PARAMETER:  7.6044E-02 -1.2255E+00 -1.2181E-01  4.3661E-01 -3.8516E-01  1.5966E-01  9.9182E-01 -3.3373E+00 -1.3939E-01 -1.5483E-01
             4.8720E-01
 GRADIENT:   1.8791E+02  2.0553E+01  3.8995E+00  3.2688E+02  2.3159E+01  3.8374E+01  8.2250E+00 -2.8311E-03  6.0495E+00  6.7645E-01
             4.8436E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1936.58767045027        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1888             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7634E-01  2.6518E-01  8.0034E-01  1.4000E+00  6.1559E-01  1.0615E+00  2.4373E+00  3.7888E-02  7.8710E-01  7.7469E-01
             1.4729E+00
 PARAMETER:  7.6056E-02 -1.2274E+00 -1.2272E-01  4.3650E-01 -3.8518E-01  1.5967E-01  9.9090E-01 -3.1731E+00 -1.3941E-01 -1.5530E-01
             4.8720E-01
 GRADIENT:   1.8793E+02  2.0222E+01  2.5451E+00  3.2628E+02  2.5090E+01  3.8375E+01  8.0763E+00 -4.2601E-03  6.0282E+00  6.5679E-01
             4.8992E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1936.59125687175        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2069
 NPARAMETR:  9.7629E-01  2.6493E-01  8.0036E-01  1.4003E+00  6.1527E-01  1.0614E+00  2.4386E+00  4.5704E-02  7.8725E-01  7.7447E-01
             1.4729E+00
 PARAMETER:  7.6009E-02 -1.2283E+00 -1.2269E-01  4.3670E-01 -3.8569E-01  1.5961E-01  9.9144E-01 -2.9856E+00 -1.3921E-01 -1.5558E-01
             4.8721E-01
 GRADIENT:   5.6755E-01  1.0609E-01  9.7116E-02 -3.3532E+00  3.0048E-01  1.5283E-01  9.3645E-02 -2.0724E-02  5.1505E-02  4.4481E-02
             3.4518E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1936.97676479087        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2253
 NPARAMETR:  9.7428E-01  2.5066E-01  8.0883E-01  1.4174E+00  6.1142E-01  1.0586E+00  2.4653E+00  4.2904E-01  7.8677E-01  7.4023E-01
             1.4626E+00
 PARAMETER:  7.3945E-02 -1.2837E+00 -1.1217E-01  4.4883E-01 -3.9197E-01  1.5698E-01  1.0023E+00 -7.4621E-01 -1.3982E-01 -2.0079E-01
             4.8023E-01
 GRADIENT:  -2.6923E+00  2.2080E+00  3.1866E-01  1.1566E+01 -9.1811E+00 -8.3073E-01 -6.9967E-01  2.5973E-01 -5.3564E-01  3.2770E+00
            -9.2898E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1937.15150843662        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2431
 NPARAMETR:  9.7496E-01  2.2272E-01  8.2740E-01  1.4295E+00  6.1610E-01  1.0598E+00  2.6962E+00  4.4109E-01  7.7677E-01  7.2176E-01
             1.4707E+00
 PARAMETER:  7.4639E-02 -1.4018E+00 -8.9461E-02  4.5735E-01 -3.8435E-01  1.5808E-01  1.0918E+00 -7.1850E-01 -1.5261E-01 -2.2607E-01
             4.8576E-01
 GRADIENT:   2.6843E-01  9.9627E-01  2.2432E+00 -1.1917E+00 -4.5905E+00 -9.9577E-02  7.9048E-02 -4.9722E-01 -6.6793E-01 -1.0023E+00
             1.4804E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1937.23283900412        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2607
 NPARAMETR:  9.7226E-01  1.4620E-01  8.8657E-01  1.4861E+00  6.2814E-01  1.0576E+00  3.4210E+00  5.1371E-01  7.5547E-01  7.3947E-01
             1.4716E+00
 PARAMETER:  7.1867E-02 -1.8228E+00 -2.0396E-02  4.9613E-01 -3.6499E-01  1.5596E-01  1.3299E+00 -5.6609E-01 -1.8042E-01 -2.0182E-01
             4.8637E-01
 GRADIENT:  -8.3804E-01  1.7409E+00  7.0855E+00  1.2639E+01 -9.8246E+00 -3.1653E-01 -1.1313E-01 -7.9271E-01 -1.2650E+00 -6.0980E-01
             1.3280E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1937.33734930846        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2786
 NPARAMETR:  9.7006E-01  9.2489E-02  9.4489E-01  1.5230E+00  6.4805E-01  1.0564E+00  4.2966E+00  6.4905E-01  7.4446E-01  7.3901E-01
             1.4630E+00
 PARAMETER:  6.9600E-02 -2.2807E+00  4.3310E-02  5.2069E-01 -3.3378E-01  1.5484E-01  1.5578E+00 -3.3224E-01 -1.9510E-01 -2.0244E-01
             4.8050E-01
 GRADIENT:  -2.1526E+00  9.7552E-01 -1.6090E+00  1.0104E+01  4.7540E+00 -4.3658E-01  2.1633E-01  3.7030E-01 -2.5619E-01 -2.5910E-01
            -1.9042E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1937.65139559699        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2963
 NPARAMETR:  9.7123E-01  5.3625E-02  9.8078E-01  1.5373E+00  6.5624E-01  1.0573E+00  5.5009E+00  6.5566E-01  7.3245E-01  7.5337E-01
             1.4708E+00
 PARAMETER:  7.0811E-02 -2.8257E+00  8.0589E-02  5.3006E-01 -3.2123E-01  1.5573E-01  1.8049E+00 -3.2212E-01 -2.1136E-01 -1.8320E-01
             4.8578E-01
 GRADIENT:   1.9819E+00 -1.6316E+00  7.0471E+00 -7.6957E+00 -1.2078E+00  3.2948E-01 -4.2265E+00  5.7606E-02  2.7891E+00  9.1183E-01
             2.7678E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1938.18180397030        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3142
 NPARAMETR:  9.6873E-01  3.8312E-02  9.1692E-01  1.5467E+00  6.2612E-01  1.0558E+00  6.4451E+00  5.6600E-01  7.3066E-01  7.4068E-01
             1.4599E+00
 PARAMETER:  6.8230E-02 -3.1620E+00  1.3270E-02  5.3615E-01 -3.6822E-01  1.5433E-01  1.9633E+00 -4.6916E-01 -2.1380E-01 -2.0019E-01
             4.7838E-01
 GRADIENT:  -1.8116E+00 -1.0418E+00  1.4939E+00  2.7923E+01 -1.7855E+00 -2.1830E-01 -2.7175E+00 -7.2743E-01  1.4961E+00 -2.2724E-02
            -4.8203E+00

0ITERATION NO.:   99    OBJECTIVE VALUE:  -1938.52965128967        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     3280
 NPARAMETR:  9.6949E-01  3.2884E-02  8.5968E-01  1.5229E+00  5.9729E-01  1.0561E+00  6.8051E+00  5.7115E-01  7.3162E-01  7.0877E-01
             1.4707E+00
 PARAMETER:  6.9004E-02 -3.3133E+00 -5.1084E-02  5.2036E-01 -4.1547E-01  1.5458E-01  2.0186E+00 -4.5554E-01 -2.1182E-01 -2.4363E-01
             4.8611E-01
 GRADIENT:  -1.4907E-01  5.8779E+01  7.1521E-01 -2.5410E+02 -3.2522E+00 -2.2402E-02  9.3602E+01  4.6805E-01  1.9247E+00  4.9232E-01
             1.8199E+00
 NUMSIGDIG:         3.1         2.3         1.9         2.3         2.5         3.4         2.3         1.0         1.5         1.6
                    2.1

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3280
 NO. OF SIG. DIGITS IN FINAL EST.:  1.0
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1038E-04  2.0276E-02 -8.5046E-03 -8.0540E-03 -1.0323E-02
 SE:             2.9758E-02  9.9317E-03  1.2353E-02  2.8216E-02  2.1057E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8095E-01  4.1195E-02  4.9115E-01  7.7530E-01  6.2398E-01

 ETASHRINKSD(%)  3.0702E-01  6.6728E+01  5.8617E+01  5.4744E+00  2.9455E+01
 ETASHRINKVR(%)  6.1310E-01  8.8929E+01  8.2875E+01  1.0649E+01  5.0234E+01
 EBVSHRINKSD(%)  6.3515E-01  7.7984E+01  5.8536E+01  4.4640E+00  2.8287E+01
 EBVSHRINKVR(%)  1.2663E+00  9.5153E+01  8.2807E+01  8.7286E+00  4.8573E+01
 RELATIVEINF(%)  9.8529E+01  2.9335E+00  1.7873E+00  4.9437E+01  5.4129E+00
 EPSSHRINKSD(%)  3.0939E+01
 EPSSHRINKVR(%)  5.2306E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1938.5296512896743     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1019.5911180850017     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.14
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1938.530       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.69E-01  3.29E-02  8.60E-01  1.52E+00  5.97E-01  1.06E+00  6.81E+00  5.74E-01  7.32E-01  7.09E-01  1.47E+00
 


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
+        1.05E+03
 
 TH 2
+       -1.17E+02  4.23E+05
 
 TH 3
+       -1.59E+01  4.58E+02  8.49E+02
 
 TH 4
+       -4.91E+00 -2.08E+04  2.08E+04  6.04E+03
 
 TH 5
+        3.56E+01 -1.58E+03 -1.65E+03 -1.68E+04  5.67E+04
 
 TH 6
+        1.32E+00  2.05E+01  1.57E+00 -5.66E+00 -3.93E+00  1.73E+02
 
 TH 7
+       -6.81E-01  2.15E+03  3.72E+00  8.39E+00 -1.18E+01  2.00E-01  2.75E+01
 
 TH 8
+       -1.28E+00 -1.56E+02 -2.44E+01  6.91E+03 -2.98E+01  1.94E-01 -9.76E-01  3.02E+01
 
 TH 9
+       -8.11E+00  2.11E+03  1.12E+02 -1.99E+02 -4.81E+04  2.24E+00  1.76E+01  2.05E+01  7.13E+02
 
 TH10
+       -3.82E+00 -3.25E+02  3.11E+01  1.05E+04 -2.48E+02  7.15E-01 -1.98E+00  4.83E+01  6.43E+01  1.60E+02
 
 TH11
+       -1.13E+01 -5.38E+01  4.86E+00  2.54E+03 -5.44E+01  2.89E+00 -1.40E-01  1.88E+01  2.93E+01  4.18E+01  2.00E+02
 
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
 #CPUT: Total CPU Time in Seconds,       63.448
Stop Time:
Wed Sep 29 21:50:55 CDT 2021
