Sun Oct 24 04:10:04 CDT 2021
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
$DATA ../../../../data/SD4/TD2/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1694.24820231650        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7333E+02 -2.3381E+01 -2.4776E+00 -3.5992E-01  5.1664E+01  4.9640E+01  5.4977E+00 -5.3594E+00  7.7408E+00  2.1773E+00
             1.2911E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.55268792444        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0336E+00  1.0160E+00  9.2760E-01  1.0474E+00  9.3250E-01  1.0055E+00  9.6007E-01  1.0539E+00  9.7159E-01  9.4682E-01
             9.6564E-01
 PARAMETER:  1.3307E-01  1.1591E-01  2.4845E-02  1.4628E-01  3.0111E-02  1.0547E-01  5.9255E-02  1.5246E-01  7.1181E-02  4.5356E-02
             6.5040E-02
 GRADIENT:   5.8092E-01  1.9059E+00 -6.0256E+00  1.2972E+01  8.8658E+00  3.7892E-01  8.5781E-01  7.1487E-01 -9.0743E-01  4.5170E+00
            -2.0717E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.86925008081        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0366E+00  1.1542E+00  8.2656E-01  9.5649E-01  9.2797E-01  1.0174E+00  9.3241E-01  1.0145E+00  1.0394E+00  8.7614E-01
             9.7259E-01
 PARAMETER:  1.3599E-01  2.4338E-01 -9.0486E-02  5.5520E-02  2.5243E-02  1.1726E-01  3.0014E-02  1.1442E-01  1.3862E-01 -3.2226E-02
             7.2203E-02
 GRADIENT:   4.4583E+00  1.4434E+01  6.7017E+00  1.0082E+01 -9.9463E+00  4.3122E+00  3.4313E+00  2.4216E-01  1.6604E+00 -1.3472E+00
             3.0948E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.93671170693        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0360E+00  1.2943E+00  5.3486E-01  8.3957E-01  8.3163E-01  1.0081E+00  8.6778E-01  4.7761E-01  1.0775E+00  7.9074E-01
             9.6850E-01
 PARAMETER:  1.3539E-01  3.5794E-01 -5.2576E-01 -7.4869E-02 -8.4373E-02  1.0811E-01 -4.1815E-02 -6.3897E-01  1.7465E-01 -1.3479E-01
             6.7995E-02
 GRADIENT:  -1.0385E+00  1.1783E+01 -4.7880E-01  8.5299E+00 -1.0772E+01 -2.4353E-01  7.0961E-01  1.2918E+00 -1.8411E+00  2.3620E+00
            -3.1628E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1701.64357039336        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0362E+00  1.5304E+00  4.5313E-01  6.9040E-01  9.1738E-01  1.0102E+00  7.5186E-01  2.4483E-01  1.2673E+00  8.1652E-01
             9.6949E-01
 PARAMETER:  1.3556E-01  5.2550E-01 -6.9157E-01 -2.7048E-01  1.3766E-02  1.1019E-01 -1.8520E-01 -1.3072E+00  3.3687E-01 -1.0271E-01
             6.9013E-02
 GRADIENT:  -4.5923E-02  1.5224E+01  1.2292E+00  1.0636E+01 -4.7661E+00  4.0754E-01 -2.8050E+00  1.8638E-01 -9.5679E-01 -2.6424E+00
            -8.8729E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.83221374304        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0358E+00  1.6877E+00  4.2037E-01  5.9179E-01  1.0008E+00  1.0084E+00  7.0558E-01  1.2668E-01  1.4522E+00  8.9272E-01
             9.7558E-01
 PARAMETER:  1.3518E-01  6.2334E-01 -7.6662E-01 -4.2461E-01  1.0075E-01  1.0841E-01 -2.4874E-01 -1.9661E+00  4.7311E-01 -1.3479E-02
             7.5274E-02
 GRADIENT:  -6.7943E-01  1.5216E+01  2.5045E+00  8.5629E+00 -6.7598E+00 -2.3513E-01 -1.6959E+00  3.2959E-02  1.1910E+00  3.6584E-01
             1.8187E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1701.97250738581        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0357E+00  1.7311E+00  3.9998E-01  5.5469E-01  1.0236E+00  1.0091E+00  7.0338E-01  4.3618E-02  1.4969E+00  8.9823E-01
             9.7010E-01
 PARAMETER:  1.3504E-01  6.4874E-01 -8.1635E-01 -4.8934E-01  1.2329E-01  1.0902E-01 -2.5186E-01 -3.0323E+00  5.0340E-01 -7.3298E-03
             6.9641E-02
 GRADIENT:  -4.9461E-01 -1.5815E+00  1.5503E+00  6.7245E-01 -1.0434E+00  4.0171E-02  1.0314E-01  4.2288E-03 -6.2617E-01 -3.3772E-03
            -4.3719E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.99597094833        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  1.0374E+00  1.7433E+00  3.8664E-01  5.4387E-01  1.0237E+00  1.0093E+00  6.9819E-01  1.0000E-02  1.5142E+00  8.9329E-01
             9.7117E-01
 PARAMETER:  1.3668E-01  6.5578E-01 -8.5027E-01 -5.0905E-01  1.2340E-01  1.0926E-01 -2.5927E-01 -5.7080E+00  5.1485E-01 -1.2845E-02
             7.0746E-02
 GRADIENT:   3.3619E+00 -7.6096E+00  1.2715E-01 -2.2902E-02  2.1943E-01  1.8666E-01  3.9829E-02  0.0000E+00  9.4947E-02  1.0213E-01
             1.6488E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1701.99646246644        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     1368
 NPARAMETR:  1.0374E+00  1.7433E+00  3.8628E-01  5.4407E-01  1.0222E+00  1.0093E+00  6.9837E-01  1.0000E-02  1.5138E+00  8.9195E-01
             9.7074E-01
 PARAMETER:  1.3676E-01  6.5524E-01 -8.5102E-01 -5.0973E-01  1.2323E-01  1.0929E-01 -2.5922E-01 -5.7080E+00  5.1489E-01 -1.3845E-02
             7.0289E-02
 GRADIENT:   5.7914E-02 -6.3584E-01  1.7096E-02 -2.8091E-01  7.7885E-01  4.0550E-03 -1.4360E-02  0.0000E+00  2.4289E-02  3.0798E-02
            -2.1733E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1368
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4502E-04 -3.5003E-02 -2.6905E-04  2.8254E-02 -3.8667E-02
 SE:             2.9859E-02  2.3396E-02  1.0154E-04  2.4247E-02  2.1971E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9345E-01  1.3464E-01  8.0584E-03  2.4391E-01  7.8414E-02

 ETASHRINKSD(%)  1.0000E-10  2.1619E+01  9.9660E+01  1.8769E+01  2.6396E+01
 ETASHRINKVR(%)  1.0000E-10  3.8564E+01  9.9999E+01  3.4015E+01  4.5825E+01
 EBVSHRINKSD(%)  4.0412E-01  2.1458E+01  9.9691E+01  1.8702E+01  2.5886E+01
 EBVSHRINKVR(%)  8.0660E-01  3.8312E+01  9.9999E+01  3.3907E+01  4.5072E+01
 RELATIVEINF(%)  9.9089E+01  4.4471E+00  9.3993E-05  5.1293E+00  9.8665E+00
 EPSSHRINKSD(%)  4.5489E+01
 EPSSHRINKVR(%)  7.0286E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.9964624664399     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.84563590270170     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.996       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.74E+00  3.86E-01  5.43E-01  1.02E+00  1.01E+00  6.98E-01  1.00E-02  1.51E+00  8.92E-01  9.71E-01
 


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
 #CPUT: Total CPU Time in Seconds,       38.917
Stop Time:
Sun Oct 24 04:10:13 CDT 2021
