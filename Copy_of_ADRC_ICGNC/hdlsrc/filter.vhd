-- -------------------------------------------------------------
--
-- Module: filter
-- Generated by MATLAB(R) 9.7 and Filter Design HDL Coder 3.1.6.
-- Generated on: 2022-01-18 23:54:20
-- -------------------------------------------------------------

-- -------------------------------------------------------------
-- HDL Code Generation Options:
--
-- TargetLanguage: VHDL
-- TestBenchStimulus: step ramp chirp 

-- Filter Specifications:
--
-- Sample Rate     : N/A (normalized frequency)
-- Response        : Lowpass
-- Specification   : N,Fp,Ap,Ast
-- Filter Order    : 2
-- Passband Ripple : 1 dB
-- Passband Edge   : 0.0008
-- Stopband Atten. : 80 dB
-- -------------------------------------------------------------

-- -------------------------------------------------------------
-- HDL Implementation    : Fully parallel
-- Folding Factor        : 1
-- -------------------------------------------------------------
-- Filter Settings:
--
-- Discrete-Time IIR Filter (real)
-- -------------------------------
-- Filter Structure    : Direct-Form II, Second-Order Sections
-- Number of Sections  : 1
-- Stable              : Yes
-- Linear Phase        : No
-- -------------------------------------------------------------



LIBRARY IEEE;
USE IEEE.std_logic_1164.all;
USE IEEE.numeric_std.ALL;

ENTITY filter IS
   PORT( clk                             :   IN    std_logic; 
         clk_enable                      :   IN    std_logic; 
         reset                           :   IN    std_logic; 
         filter_in                       :   IN    real; -- double
         filter_out                      :   OUT   real  -- double
         );

END filter;


----------------------------------------------------------------
--Module Architecture: filter
----------------------------------------------------------------
ARCHITECTURE rtl OF filter IS
  -- Local Functions
  -- Type Definitions
  TYPE delay_pipeline_type IS ARRAY (NATURAL range <>) OF real; -- double
  -- Constants
  CONSTANT scaleconst1                    : real := 1.0141169989119779E-04; -- double
  CONSTANT coeff_b1_section1              : real := 1.0000000000000000E+00; -- double
  CONSTANT coeff_b2_section1              : real := -1.9388780109630124E+00; -- double
  CONSTANT coeff_b3_section1              : real := 1.0000000000000000E+00; -- double
  CONSTANT coeff_a2_section1              : real := -1.9972380961558445E+00; -- double
  CONSTANT coeff_a3_section1              : real := 9.9724505097018890E-01; -- double
  -- Signals
  SIGNAL input_register                   : real := 0.0; -- double
  SIGNAL scale1                           : real := 0.0; -- double
  SIGNAL scaletypeconvert1                : real := 0.0; -- double
  -- Section 1 Signals 
  SIGNAL a1sum1                           : real := 0.0; -- double
  SIGNAL a2sum1                           : real := 0.0; -- double
  SIGNAL b1sum1                           : real := 0.0; -- double
  SIGNAL b2sum1                           : real := 0.0; -- double
  SIGNAL delay_section1                   : delay_pipeline_type(0 TO 1) := (0.0, 0.0); -- double
  SIGNAL inputconv1                       : real := 0.0; -- double
  SIGNAL a2mul1                           : real := 0.0; -- double
  SIGNAL a3mul1                           : real := 0.0; -- double
  SIGNAL b1mul1                           : real := 0.0; -- double
  SIGNAL b2mul1                           : real := 0.0; -- double
  SIGNAL b3mul1                           : real := 0.0; -- double
  SIGNAL output_typeconvert               : real := 0.0; -- double
  SIGNAL output_register                  : real := 0.0; -- double


BEGIN

  -- Block Statements
  input_reg_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      input_register <= 0.0000000000000000E+00;
    ELSIF clk'event AND clk = '1' THEN
      IF clk_enable = '1' THEN
        input_register <= filter_in;
      END IF;
    END IF; 
  END PROCESS input_reg_process;

  scale1 <= input_register * scaleconst1;

  scaletypeconvert1 <= scale1;


  --   ------------------ Section 1 ------------------

  delay_process_section1 : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      delay_section1(0) <= 0.0000000000000000E+00;
      delay_section1(1) <= 0.0000000000000000E+00;
    ELSIF clk'event AND clk = '1' THEN
      IF clk_enable = '1' THEN
        delay_section1(1) <= delay_section1(0);
        delay_section1(0) <= a1sum1;
      END IF;
    END IF;
  END PROCESS delay_process_section1;

  inputconv1 <= scaletypeconvert1;


  a2mul1 <= delay_section1(0) * coeff_a2_section1;

  a3mul1 <= delay_section1(1) * coeff_a3_section1;

  b1mul1 <= a1sum1;


  b2mul1 <= delay_section1(0) * coeff_b2_section1;

  b3mul1 <= delay_section1(1);


  a2sum1 <= inputconv1 - a2mul1;

  a1sum1 <= a2sum1 - a3mul1;

  b2sum1 <= b1mul1 + b2mul1;

  b1sum1 <= b2sum1 + b3mul1;

  output_typeconvert <= b1sum1;


  Output_Register_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      output_register <= 0.0000000000000000E+00;
    ELSIF clk'event AND clk = '1' THEN
      IF clk_enable = '1' THEN
        output_register <= output_typeconvert;
      END IF;
    END IF; 
  END PROCESS Output_Register_process;

  -- Assignment Statements
  filter_out <= output_register;
END rtl;
