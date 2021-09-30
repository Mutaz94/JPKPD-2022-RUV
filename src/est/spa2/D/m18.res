Thu Sep 30 08:39:57 CDT 2021
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
$DATA ../../../../data/spa2/D/dat18.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m18.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1342.42381399569        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.9203E+01 -5.4924E+01 -7.3154E+01 -1.8642E+02  4.5689E+02 -7.9495E+02 -2.8816E+02 -5.5302E+01 -7.1631E+02 -5.4393E+02
             6.2451E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2169.88384017653        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.2249E+00  9.7477E-01  8.9850E-01  1.2688E+00  8.2071E-01  2.3807E+00  1.6129E+00  1.2061E+00  3.0222E+00  2.2328E+00
             9.1763E-01
 PARAMETER:  3.0282E-01  7.4443E-02 -7.0264E-03  3.3806E-01 -9.7583E-02  9.6741E-01  5.7805E-01  2.8741E-01  1.2060E+00  9.0325E-01
             1.4044E-02
 GRADIENT:   3.5221E+01 -3.5470E+01 -3.0715E+01  1.6207E+01  2.5379E+01  3.6017E+01  6.0702E+00 -7.9707E+00  4.2564E+01  2.2882E+01
            -1.8187E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2185.14767886927        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  1.1064E+00  1.1618E+00  1.9516E+00  1.4859E+00  1.2980E+00  2.5440E+00  2.9055E+00  2.6623E+00  3.0538E+00  1.7993E+00
             9.1580E-01
 PARAMETER:  2.0107E-01  2.4996E-01  7.6867E-01  4.9601E-01  3.6080E-01  1.0337E+00  1.1666E+00  1.0792E+00  1.2164E+00  6.8739E-01
             1.2042E-02
 GRADIENT:  -9.9651E+00 -7.0542E+00 -2.6264E+01  4.7781E+01  3.9284E+01  1.6782E+02 -1.3556E+01  1.2348E+01  4.1915E+01  4.1666E+01
            -2.7774E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2267.67764668737        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.2010E+00  7.5475E-01  3.7821E+00  1.1861E+00  1.3560E+00  2.4539E+00  4.9012E+00  2.8001E+00  1.6124E+00  1.5048E+00
             9.8037E-01
 PARAMETER:  2.8317E-01 -1.8137E-01  1.4303E+00  2.7069E-01  4.0452E-01  9.9770E-01  1.6895E+00  1.1296E+00  5.7773E-01  5.0865E-01
             8.0178E-02
 GRADIENT:   1.2956E+03  7.9151E+01  5.1243E+01  3.4742E+02  3.5492E+01  2.3907E+03  3.2673E+03  4.6208E+01  1.7003E+02  2.2279E+01
             5.6204E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2273.23335192147        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  1.1630E+00  8.2659E-01  3.7934E+00  1.3195E+00  1.3624E+00  2.3697E+00  5.1958E+00  2.7382E+00  1.3257E+00  1.5150E+00
             9.7469E-01
 PARAMETER:  2.5098E-01 -9.0446E-02  1.4333E+00  3.7727E-01  4.0923E-01  9.6275E-01  1.7479E+00  1.1073E+00  3.8196E-01  5.1543E-01
             7.4362E-02
 GRADIENT:   1.0589E+01  5.2504E+00  1.9323E+01  1.0619E+01 -2.0861E+01  1.2656E+02 -1.6623E+02  1.0948E+00 -2.7576E+00  2.4542E+00
            -2.0339E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2275.26627839803        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      876             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1683E+00  8.6341E-01  2.6180E+00  1.2760E+00  1.3683E+00  2.3958E+00  5.1334E+00  2.4659E+00  1.2406E+00  1.4990E+00
             9.8001E-01
 PARAMETER:  2.5557E-01 -4.6861E-02  1.0624E+00  3.4369E-01  4.1358E-01  9.7374E-01  1.7358E+00  1.0026E+00  3.1560E-01  5.0479E-01
             7.9810E-02
 GRADIENT:   1.1537E+03  5.0787E+01  2.0739E+01  5.3253E+02  1.0657E+02  2.3522E+03  3.0385E+03  3.0363E+01  5.8548E+01  3.1362E+01
             1.0601E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2275.97205722509        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      992
 NPARAMETR:  1.1916E+00  8.4780E-01  2.5607E+00  1.2513E+00  1.3094E+00  2.2240E+00  5.0231E+00  2.3659E+00  1.2716E+00  1.4744E+00
             9.7473E-01
 PARAMETER:  2.7528E-01 -6.5115E-02  1.0403E+00  3.2422E-01  3.6957E-01  8.9930E-01  1.7140E+00  9.6118E-01  3.4031E-01  4.8824E-01
             7.4401E-02
 GRADIENT:   2.3505E+01 -1.5593E+00  1.2695E+00  1.8478E+01  4.7260E+00  7.7248E+01 -1.5360E+02  1.1936E+00  3.8284E+00  5.7479E+00
             5.4449E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2278.31584773288        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1133
 NPARAMETR:  1.1578E+00  8.9584E-01  2.4296E+00  1.1650E+00  1.2971E+00  2.3420E+00  4.6815E+00  2.3472E+00  1.2194E+00  1.4066E+00
             9.6836E-01
 PARAMETER:  2.4652E-01 -9.9910E-03  9.8773E-01  2.5268E-01  3.6013E-01  9.5103E-01  1.6436E+00  9.5322E-01  2.9837E-01  4.4116E-01
             6.7851E-02
 GRADIENT:   8.3375E+00 -7.5209E+00  7.7176E+00 -4.3363E+00 -1.0979E+01  1.2142E+02 -1.4759E+02  2.7288E+00  2.9341E+00 -4.8719E+00
            -6.3430E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2279.09315023680        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.1629E+00  9.9294E-01  2.2907E+00  1.1661E+00  1.3153E+00  2.3380E+00  4.7049E+00  2.3076E+00  1.1867E+00  1.4411E+00
             9.7525E-01
 PARAMETER:  2.5093E-01  9.2920E-02  9.2885E-01  2.5368E-01  3.7410E-01  9.4929E-01  1.6486E+00  9.3623E-01  2.7119E-01  4.6544E-01
             7.4937E-02
 GRADIENT:   1.1421E+03  5.7250E+01  2.3436E+01  3.4302E+02  7.2741E+01  2.3115E+03  2.7630E+03  2.3358E+01  4.1285E+01  2.2502E+01
             2.2109E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2279.67828611368        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  1.1697E+00  1.0226E+00  2.2380E+00  1.0929E+00  1.3072E+00  2.3438E+00  4.6699E+00  2.2288E+00  1.1294E+00  1.4290E+00
             9.7268E-01
 PARAMETER:  2.5672E-01  1.2232E-01  9.0558E-01  1.8881E-01  3.6791E-01  9.5176E-01  1.6411E+00  9.0147E-01  2.2166E-01  4.5695E-01
             7.2300E-02
 GRADIENT:   1.3154E+01 -3.4089E+00  9.3795E+00 -1.0067E+01 -1.3217E+01  1.1839E+02 -1.2122E+02 -2.1022E+00  6.6239E-01 -1.1500E+00
            -7.4181E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2280.11669270278        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1540
 NPARAMETR:  1.1578E+00  1.0667E+00  2.0884E+00  1.1001E+00  1.3267E+00  2.3356E+00  4.4204E+00  2.2498E+00  1.1344E+00  1.4308E+00
             9.7951E-01
 PARAMETER:  2.4656E-01  1.6461E-01  8.3642E-01  1.9544E-01  3.8269E-01  9.4828E-01  1.5862E+00  9.1083E-01  2.2614E-01  4.5825E-01
             7.9297E-02
 GRADIENT:   7.9050E+00 -2.7503E+00 -3.6485E+00  1.8323E+01  6.7686E+00  1.1673E+02 -1.1656E+02  3.3388E+00 -5.2496E-01 -3.8052E-01
             5.5599E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2280.77022874611        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1730             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1646E+00  1.1268E+00  2.0648E+00  1.0505E+00  1.3230E+00  2.3768E+00  4.4150E+00  2.1999E+00  1.1193E+00  1.4344E+00
             9.8057E-01
 PARAMETER:  2.5240E-01  2.1942E-01  8.2505E-01  1.4927E-01  3.7992E-01  9.6575E-01  1.5850E+00  8.8842E-01  2.1272E-01  4.6072E-01
             8.0378E-02
 GRADIENT:   1.1379E+03  1.4185E+02  2.3825E+01  1.6051E+02  6.3110E+01  2.3363E+03  2.5383E+03  1.5980E+01  2.4922E+01  2.0625E+01
            -6.5352E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2280.85566760248        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:     1854
 NPARAMETR:  1.1540E+00  1.1802E+00  1.9346E+00  1.0403E+00  1.3431E+00  2.3734E+00  4.2870E+00  2.2031E+00  1.1199E+00  1.4336E+00
             9.8292E-01
 PARAMETER:  2.4321E-01  2.6565E-01  7.5993E-01  1.3949E-01  3.9497E-01  9.6433E-01  1.5556E+00  8.8987E-01  2.1324E-01  4.6022E-01
             8.2772E-02
 GRADIENT:   1.0865E+03  1.7817E+02  1.2723E+01  1.5945E+02  8.2666E+01  2.3337E+03  2.4109E+03  1.7801E+01  2.4419E+01  2.0721E+01
             1.9712E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2281.07462156986        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1999             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1746E+00  1.2077E+00  1.9538E+00  1.0050E+00  1.3321E+00  2.4093E+00  4.2347E+00  2.1476E+00  1.0769E+00  1.4380E+00
             9.8330E-01
 PARAMETER:  2.6092E-01  2.8869E-01  7.6980E-01  1.0500E-01  3.8678E-01  9.7935E-01  1.5433E+00  8.6433E-01  1.7408E-01  4.6326E-01
             8.3156E-02
 GRADIENT:   1.1735E+03  1.9937E+02  2.2932E+01  1.0166E+02  6.3436E+01  2.3523E+03  2.3844E+03  1.2106E+01  1.5249E+01  2.1264E+01
            -1.1859E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2281.18437773892        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2161
 NPARAMETR:  1.1659E+00  1.2231E+00  1.8886E+00  1.0113E+00  1.3405E+00  2.3760E+00  4.2430E+00  2.1772E+00  1.0932E+00  1.4363E+00
             9.8427E-01
 PARAMETER:  2.5348E-01  3.0135E-01  7.3582E-01  1.1126E-01  3.9301E-01  9.6543E-01  1.5453E+00  8.7804E-01  1.8911E-01  4.6208E-01
             8.4149E-02
 GRADIENT:   1.1048E+01 -5.7720E-01 -2.5917E+00  7.2670E+00  4.9160E+00  1.2549E+02 -9.3172E+01  1.8361E+00  2.9296E-01  2.0834E-01
             5.5648E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2281.26895370159        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2354             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1653E+00  1.2364E+00  1.8905E+00  9.9529E-01  1.3350E+00  2.3781E+00  4.2097E+00  2.1523E+00  1.0902E+00  1.4354E+00
             9.8421E-01
 PARAMETER:  2.5298E-01  3.1217E-01  7.3686E-01  9.5278E-02  3.8890E-01  9.6628E-01  1.5374E+00  8.6653E-01  1.8639E-01  4.6145E-01
             8.4085E-02
 GRADIENT:   1.1322E+03  2.2021E+02  1.8515E+01  9.6722E+01  6.7574E+01  2.3204E+03  2.3430E+03  1.3709E+01  1.8043E+01  2.0806E+01
             7.4078E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2281.28801392249        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2546
 NPARAMETR:  1.1653E+00  1.2467E+00  1.8741E+00  9.8759E-01  1.3360E+00  2.3783E+00  4.1845E+00  2.1460E+00  1.0891E+00  1.4350E+00
             9.8454E-01
 PARAMETER:  2.5294E-01  3.2046E-01  7.2812E-01  8.7512E-02  3.8967E-01  9.6640E-01  1.5314E+00  8.6362E-01  1.8533E-01  4.6115E-01
             8.4417E-02
 GRADIENT:   1.0776E+01 -1.4371E+00  1.8645E+00 -2.7899E+00 -2.9175E+00  1.2626E+02 -9.0520E+01 -9.8731E-02  3.2902E-01  9.1526E-02
            -5.4567E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2281.31497951868        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2739             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1653E+00  1.2672E+00  1.8436E+00  9.8710E-01  1.3393E+00  2.3783E+00  4.1563E+00  2.1431E+00  1.0852E+00  1.4353E+00
             9.8511E-01
 PARAMETER:  2.5298E-01  3.3679E-01  7.1170E-01  8.7021E-02  3.9212E-01  9.6639E-01  1.5246E+00  8.6227E-01  1.8173E-01  4.6138E-01
             8.5000E-02
 GRADIENT:   1.1299E+03  2.4211E+02  1.5717E+01  9.8817E+01  7.1361E+01  2.3163E+03  2.2857E+03  1.3637E+01  1.6697E+01  2.0893E+01
             1.2371E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2281.32387638507        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2931
 NPARAMETR:  1.1653E+00  1.2728E+00  1.8357E+00  9.8325E-01  1.3398E+00  2.3783E+00  4.1435E+00  2.1391E+00  1.0846E+00  1.4353E+00
             9.8530E-01
 PARAMETER:  2.5300E-01  3.4122E-01  7.0745E-01  8.3107E-02  3.9249E-01  9.6639E-01  1.5216E+00  8.6041E-01  1.8122E-01  4.6139E-01
             8.5195E-02
 GRADIENT:   1.0756E+01 -2.5175E-01 -3.1663E-01  2.1795E+00  4.5997E-01  1.2605E+02 -8.8400E+01  3.8256E-01 -1.6330E-01  1.8654E-01
             2.3846E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2281.33134421496        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3124             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1653E+00  1.2763E+00  1.8311E+00  9.7791E-01  1.3397E+00  2.3783E+00  4.1341E+00  2.1350E+00  1.0858E+00  1.4353E+00
             9.8545E-01
 PARAMETER:  2.5302E-01  3.4395E-01  7.0494E-01  7.7664E-02  3.9242E-01  9.6639E-01  1.5193E+00  8.5845E-01  1.8227E-01  4.6137E-01
             8.5339E-02
 GRADIENT:   1.1295E+03  2.4959E+02  1.6226E+01  9.4450E+01  6.9975E+01  2.3149E+03  2.2736E+03  1.3051E+01  1.6715E+01  2.0766E+01
             1.0934E+00

0ITERATION NO.:   97    OBJECTIVE VALUE:  -2281.33134421496        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     3187
 NPARAMETR:  1.1653E+00  1.2763E+00  1.8311E+00  9.7791E-01  1.3397E+00  2.3783E+00  4.1341E+00  2.1350E+00  1.0858E+00  1.4353E+00
             9.8545E-01
 PARAMETER:  2.5302E-01  3.4395E-01  7.0494E-01  7.7664E-02  3.9242E-01  9.6639E-01  1.5193E+00  8.5845E-01  1.8227E-01  4.6137E-01
             8.5339E-02
 GRADIENT:  -2.3986E-03 -1.1521E-01  3.8872E-01 -1.6622E-02 -3.2136E-01  3.5702E-04  3.4763E-01  4.3952E-02  5.1698E-02 -1.2014E-02
            -1.0601E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3187
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7834E-04  1.1786E-02 -5.5412E-02 -1.8859E-02 -2.5640E-02
 SE:             2.9967E-02  2.6801E-02  1.6409E-02  1.8937E-02  2.4240E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9259E-01  6.6012E-01  7.3317E-04  3.1929E-01  2.9017E-01

 ETASHRINKSD(%)  1.0000E-10  1.0213E+01  4.5028E+01  3.6559E+01  1.8793E+01
 ETASHRINKVR(%)  1.0000E-10  1.9383E+01  6.9781E+01  5.9753E+01  3.4054E+01
 EBVSHRINKSD(%)  6.1708E-02  6.9232E+00  4.8248E+01  4.3607E+01  1.5663E+01
 EBVSHRINKVR(%)  1.2338E-01  1.3367E+01  7.3217E+01  6.8199E+01  2.8873E+01
 RELATIVEINF(%)  9.9873E+01  5.2567E+01  1.3996E+01  1.3562E+01  4.8249E+01
 EPSSHRINKSD(%)  3.1412E+01
 EPSSHRINKVR(%)  5.2957E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2281.3313442149624     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1178.6051043693553     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    78.39
 Elapsed covariance  time in seconds:    11.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2281.331       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.17E+00  1.28E+00  1.83E+00  9.78E-01  1.34E+00  2.38E+00  4.13E+00  2.13E+00  1.09E+00  1.44E+00  9.85E-01
 


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
 
         8.31E-02  2.74E-01  4.03E-01  1.47E-01  1.00E-01  2.63E-01  5.61E-01  2.47E-01  1.80E-01  1.26E-01  4.59E-02
 


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
+        6.91E-03
 
 TH 2
+        6.30E-03  7.50E-02
 
 TH 3
+        6.14E-03 -6.89E-02  1.63E-01
 
 TH 4
+        2.45E-03 -3.12E-02  4.60E-02  2.15E-02
 
 TH 5
+        8.26E-04  1.70E-02 -3.70E-03 -7.97E-03  1.01E-02
 
 TH 6
+        1.30E-02  3.15E-03  2.96E-02  5.91E-03  4.43E-03  6.90E-02
 
 TH 7
+        1.88E-02 -8.78E-02  1.54E-01  5.50E-02 -1.24E-02  1.02E-01  3.15E-01
 
 TH 8
+        1.97E-03  2.17E-02 -7.04E-03 -1.07E-02  1.10E-02  1.70E-02 -7.69E-03  6.11E-02
 
 TH 9
+       -1.54E-03 -1.67E-02  3.20E-02  1.26E-02 -3.72E-03 -3.89E-03  1.42E-02 -4.54E-03  3.23E-02
 
 TH10
+        1.54E-03  7.17E-04  5.39E-03  1.64E-03 -4.46E-04  3.34E-03  1.26E-02  5.20E-03  5.25E-04  1.60E-02
 
 TH11
+        1.79E-05  3.52E-03 -5.61E-03 -2.71E-03  4.94E-04 -6.67E-04 -5.82E-03  2.03E-03 -1.70E-03 -1.41E-04  2.11E-03
 
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
+        8.31E-02
 
 TH 2
+        2.77E-01  2.74E-01
 
 TH 3
+        1.83E-01 -6.24E-01  4.03E-01
 
 TH 4
+        2.01E-01 -7.78E-01  7.78E-01  1.47E-01
 
 TH 5
+        9.91E-02  6.21E-01 -9.16E-02 -5.42E-01  1.00E-01
 
 TH 6
+        5.96E-01  4.39E-02  2.80E-01  1.53E-01  1.68E-01  2.63E-01
 
 TH 7
+        4.04E-01 -5.72E-01  6.80E-01  6.69E-01 -2.20E-01  6.95E-01  5.61E-01
 
 TH 8
+        9.57E-02  3.21E-01 -7.07E-02 -2.96E-01  4.44E-01  2.62E-01 -5.55E-02  2.47E-01
 
 TH 9
+       -1.03E-01 -3.39E-01  4.42E-01  4.80E-01 -2.07E-01 -8.23E-02  1.41E-01 -1.02E-01  1.80E-01
 
 TH10
+        1.46E-01  2.07E-02  1.06E-01  8.82E-02 -3.52E-02  1.00E-01  1.77E-01  1.67E-01  2.31E-02  1.26E-01
 
 TH11
+        4.70E-03  2.80E-01 -3.03E-01 -4.03E-01  1.07E-01 -5.53E-02 -2.26E-01  1.79E-01 -2.06E-01 -2.42E-02  4.59E-02
 
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
+        6.20E+02
 
 TH 2
+       -2.39E+02  1.70E+02
 
 TH 3
+       -1.39E+01  2.12E+01  3.38E+01
 
 TH 4
+       -2.77E+02  9.64E+01 -5.90E+01  4.79E+02
 
 TH 5
+        1.09E+02 -1.29E+02 -7.37E+01  1.55E+02  3.88E+02
 
 TH 6
+       -1.78E+01 -3.96E+01 -8.06E-01  3.80E+01  2.53E+01  7.86E+01
 
 TH 7
+       -4.36E+01  4.50E+01 -1.82E+00 -1.16E+01 -2.74E+01 -4.12E+01  3.40E+01
 
 TH 8
+        2.21E-01  7.43E+00  5.64E-02  1.18E+00 -2.51E+01 -1.35E+01  6.01E+00  2.44E+01
 
 TH 9
+        5.18E+01 -1.93E+01 -8.74E+00 -5.40E+01  1.18E+01 -5.71E+00  5.04E+00 -1.14E+00  5.11E+01
 
 TH10
+        2.25E+01 -3.39E+01 -5.32E+00  1.87E+00  4.00E+01  2.14E+01 -1.71E+01 -1.12E+01  5.38E-01  7.73E+01
 
 TH11
+       -1.09E+02  1.64E+01 -1.66E+01  2.00E+02  9.46E+01  2.79E+01 -9.94E+00 -1.77E+01 -9.53E+00  1.11E+01  6.31E+02
 
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
 #CPUT: Total CPU Time in Seconds,       89.483
Stop Time:
Thu Sep 30 08:41:29 CDT 2021
