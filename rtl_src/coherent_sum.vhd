library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use IEEE.math_real.all;
use work.radar_pkg.all;

entity coherent_sum is
    Port (
        clk_i   : in std_logic;
        rst_i   : in std_logic;
        data_i  : in std_logic_vector (ADC_BIT_RES_C-1 downto 0);
        ena_i   : in std_logic;
        sum_o   : out std_logic_vector (ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C-1 downto 0);
        vld_o   : out std_logic
    );
end coherent_sum;

architecture rtl of coherent_sum is

    attribute ram_style : string;

    constant coh_pls_c  : integer := integer(COHERENT_NP_C);
    type acc_ram_t is array(0 to PULSE_RXSAMPLES_C-1) of signed(ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C-1 downto 0);
    signal acc_ram_s    : acc_ram_t := (others => (others => '0'));
    attribute ram_style of acc_ram_s : signal is "block";

    signal rd_ctr_s     : integer range 0 to coh_pls_c - 1 := 0;
    signal wr_ctr_s     : integer range 0 to coh_pls_c - 1 := 0;
    signal ext_bin_s    : signed(ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C-1 downto 0) := (others => '0');
    signal rdata_s      : signed(ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C-1 downto 0);
    signal slow_ctr_s   : integer range 0 to integer(COHERENT_NP_C) - 1 := 0;    -- slow domain counter up to Np!

    signal vld_s        : std_logic;
    signal ena_s        : std_logic_vector(1 downto 0);
    signal wdata_s      : signed(ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C-1 downto 0):= (others => '0');

begin

    fdbk_acc: process(clk_i)
        variable edge_v     : std_logic:= '0';
    begin
        if rising_edge(clk_i) then
            if rst_i = '1' then
                wr_ctr_s   <= 0;
                rd_ctr_s    <= 0;
                slow_ctr_s  <= 0;
                vld_s       <= '0';
                ena_s       <= (others => '0');
                ext_bin_s   <= (others => '0');
            else

                vld_s       <= '0';
                rdata_s     <= acc_ram_s(rd_ctr_s);
                ena_s       <= ena_s(0) & ena_i;
                ext_bin_s   <= resize(signed(data_i), ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C);
                if ena_i = '1' then                         -- read from RAM
                    if rd_ctr_s < PULSE_RXSAMPLES_C - 1 then
                        rd_ctr_s       <= rd_ctr_s + 1;     -- generate rd_ctr before data mux for write
                    end if;
                elsif ena_s(0) = '0' and ena_s(1) = '1' then
                    if slow_ctr_s < coh_pls_c - 1 then
                        slow_ctr_s  <= slow_ctr_s + 1;
                    else
                        slow_ctr_s  <= 0;
                    end if;
                    rd_ctr_s       <= 0;
                end if;
                if ena_s(0) = '1' then                      -- provide data for write
                    if slow_ctr_s = 0 then      -- load each bin in the respective cell
                        wdata_s     <= ext_bin_s;
                    elsif slow_ctr_s <= coh_pls_c - 1 then
                        wdata_s     <= rdata_s + ext_bin_s;
                        if slow_ctr_s = coh_pls_c - 1 then
                            vld_s   <= '1';
                        end if;
                    end if;
                end if;
                if ena_s(1) = '1' then
                    acc_ram_s(wr_ctr_s) <= wdata_s;
                    wr_ctr_s            <= wr_ctr_s + 1;
                else
                    wr_ctr_s            <= 0;
                end if;
            end if;
        end if;
    end process;

    -- register outputs
    reg_out: process(clk_i)
    begin
        if rising_edge(clk_i) then
            if rst_i = '1' then
                sum_o       <= (others => '0');
                vld_o       <= '0';
            else
                if vld_s = '1' then
                    vld_o   <= '1';
                    sum_o   <= std_logic_vector(wdata_s);
                else
                    sum_o   <= (others => '0');
                    vld_o   <= '0';
                end if;
            end if;
        end if;
    end process;

end rtl;
