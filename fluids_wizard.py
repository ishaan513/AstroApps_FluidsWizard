import streamlit as st
import numpy as np
import pandas as pd

# Fluid Properties Database
FLUIDS = {
    "LOX (Liquid Oxygen)": {
        "type": "liquid",
        "density": 1141,  # kg/m³ at NBP
        "temp": -183,  # °C
        "molecular_weight": 32.0,
        "gamma": 1.4
    },
    "RP-1 (Kerosene)": {
        "type": "liquid",
        "density": 820,  # kg/m³
        "temp": 15,  # °C
        "molecular_weight": 170,
        "gamma": None
    },
    "Liquid Methane (LCH4)": {
        "type": "liquid",
        "density": 422,  # kg/m³ at NBP
        "temp": -161,  # °C
        "molecular_weight": 16.04,
        "gamma": 1.31
    },
    "Ethanol": {
        "type": "liquid",
        "density": 789,  # kg/m³
        "temp": 20,  # °C
        "molecular_weight": 46.07,
        "gamma": None
    },
    "Water": {
        "type": "liquid",
        "density": 1000,  # kg/m³
        "temp": 20,  # °C
        "molecular_weight": 18.015,
        "gamma": None
    },
    "Nitrous Oxide (N2O)": {
        "type": "liquid",
        "density": 1220,  # kg/m³ at 20°C
        "temp": 20,  # °C
        "molecular_weight": 44.013,
        "gamma": 1.27
    },
    "Isopropyl Alcohol (IPA)": {
        "type": "liquid",
        "density": 786,  # kg/m³ at 20°C
        "temp": 20,  # °C
        "molecular_weight": 60.1,
        "gamma": None
    },
    "Helium": {
        "type": "gas",
        "density": None,  # Calculate with ideal gas
        "temp": 15,  # °C (reference)
        "molecular_weight": 4.003,
        "gamma": 1.66,
        "R_specific": 2077  # J/(kg·K)
    },
    "Nitrogen": {
        "type": "gas",
        "density": None,
        "temp": 15,  # °C
        "molecular_weight": 28.014,
        "gamma": 1.4,
        "R_specific": 297  # J/(kg·K)
    },
    "Gaseous Oxygen": {
        "type": "gas",
        "density": None,
        "temp": 15,  # °C
        "molecular_weight": 32.0,
        "gamma": 1.4,
        "R_specific": 260  # J/(kg·K)
    },
    "Gaseous Methane (GCH4)": {
        "type": "gas",
        "density": None,
        "temp": 15,  # °C
        "molecular_weight": 16.04,
        "gamma": 1.31,
        "R_specific": 518  # J/(kg·K)
    }
}

st.set_page_config(page_title="Astro Fluids Wizard", page_icon="🚀", layout="wide")

st.title("🚀 Astro Fluids Wizard")
st.markdown("*Fluid dynamics tool*")

# Sidebar for fluid properties reference
with st.sidebar:
    st.header("Fluid Properties Reference")
    fluid_ref = st.selectbox("Select fluid:", list(FLUIDS.keys()))
    props = FLUIDS[fluid_ref]
    st.write(f"**Type:** {props['type']}")
    if props['density']:
        st.write(f"**Density:** {props['density']} kg/m³")
    st.write(f"**Temp (ref):** {props['temp']} °C")
    st.write(f"**Molecular Weight:** {props['molecular_weight']} g/mol")
    if props['gamma']:
        st.write(f"**Specific Heat Ratio (γ):** {props['gamma']}")
    if props['type'] == 'gas':
        st.write(f"**R (specific):** {props['R_specific']} J/(kg·K)")
    
    st.markdown("---")
    
    st.header("CdA ↔ Cv Converter")
    st.markdown("*Flow coefficient conversion*")
    
    conversion_type = st.radio("Convert:", ["CdA to Cv", "Cv to CdA"], key="conv_type")
    
    if conversion_type == "CdA to Cv":
        cda_input = st.number_input("CdA (in²):", value=0.01, min_value=0.0, step=0.001, format="%.4f", key="cda_in")
        # Cv = CdA * 38.0 (approximate conversion factor for water at 60°F)
        cv_result = cda_input * 38.0
        st.success(f"**Cv = {cv_result:.3f}**")
    else:
        cv_input = st.number_input("Cv:", value=0.38, min_value=0.0, step=0.01, format="%.3f", key="cv_in")
        # CdA = Cv / 38.0
        cda_result = cv_input / 38.0
        st.success(f"**CdA = {cda_result:.4f} in²**")
    
    st.caption("*Cv: gpm at 1 psi drop*")

# Main tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "💧 Incompressible Flow", 
    "💨 Compressible Flow", 
    "⚗️ Ideal Gas Calculator",
    "🎯 Injector Sizing",
    "📊 Mixture Ratio"
])

# Tab 1: Incompressible Flow
with tab1:
    st.header("Incompressible Flow (Liquids)")
    st.markdown("Calculate mass flow rate through orifices and injectors")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Input Parameters")
        fluid_incomp = st.selectbox("Select fluid:", 
                                     [f for f in FLUIDS.keys() if FLUIDS[f]['type'] == 'liquid'],
                                     key="incomp_fluid")
        
        use_custom_density = st.checkbox("Use custom density", key="custom_dens_incomp")
        if use_custom_density:
            density_incomp = st.number_input("Density (kg/m³):", value=1000.0, min_value=0.1)
        else:
            density_incomp = FLUIDS[fluid_incomp]['density']
            st.info(f"Using density: {density_incomp} kg/m³")
        
        cd = st.number_input("Discharge coefficient (Cd):", value=0.7, min_value=0.1, max_value=1.0, step=0.05)
        diameter = st.number_input("Orifice diameter (inches):", value=0.079, min_value=0.001, step=0.001, format="%.3f")
        delta_p = st.number_input("Pressure drop (psi):", value=145.0, min_value=0.0, step=5.0)
    
    with col2:
        st.subheader("Results")
        
        # Calculate area in m²
        area = np.pi * (diameter * 0.0254 / 2)**2  # Convert inches to meters
        
        # Calculate velocity: v = sqrt(2*ΔP/ρ)
        delta_p_pa = delta_p * 6894.76  # Convert psi to Pa
        velocity = np.sqrt(2 * delta_p_pa / density_incomp)
        
        # Mass flow rate: ṁ = Cd * A * ρ * v
        mass_flow = cd * area * density_incomp * velocity
        
        # Volume flow rate
        vol_flow = mass_flow / density_incomp * 1000  # L/s
        
        st.metric("Mass Flow Rate", f"{mass_flow:.4f} kg/s")
        st.metric("Mass Flow Rate", f"{mass_flow*1000:.2f} g/s")
        st.metric("Volume Flow Rate", f"{vol_flow:.2f} L/s")
        st.metric("Velocity", f"{velocity:.2f} m/s")
        st.metric("Orifice Area", f"{area*1550:.3f} in²")

# Tab 2: Compressible Flow
with tab2:
    st.header("Compressible Flow (Gases)")
    st.markdown("Calculate mass flow rate with choked and unchoked conditions")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Input Parameters")
        fluid_comp = st.selectbox("Select fluid:", 
                                   [f for f in FLUIDS.keys() if FLUIDS[f]['type'] == 'gas'],
                                   key="comp_fluid")
        
        gamma = FLUIDS[fluid_comp]['gamma']
        R_specific = FLUIDS[fluid_comp]['R_specific']
        
        st.info(f"γ = {gamma}, R = {R_specific} J/(kg·K)")
        
        cd_comp = st.number_input("Discharge coefficient (Cd):", value=0.7, min_value=0.1, max_value=1.0, step=0.05, key="cd_comp")
        diameter_comp = st.number_input("Orifice diameter (inches):", value=0.079, min_value=0.001, step=0.001, format="%.3f", key="diam_comp")
        p_upstream = st.number_input("Upstream pressure (psi):", value=725.0, min_value=0.1, step=10.0)
        p_downstream = st.number_input("Downstream pressure (psi):", value=145.0, min_value=0.0, step=10.0)
        temp_upstream = st.number_input("Upstream temperature (°C):", value=20.0, step=1.0)
    
    with col2:
        st.subheader("Results")
        
        # Convert units
        p1 = p_upstream * 6894.76  # psi to Pa
        p2 = p_downstream * 6894.76  # psi to Pa
        T1 = temp_upstream + 273.15  # K
        area_comp = np.pi * (diameter_comp * 0.0254 / 2)**2  # inches to m²
        
        # Calculate upstream density
        rho1 = p1 / (R_specific * T1)
        
        # Critical pressure ratio
        pressure_ratio = p2 / p1
        critical_ratio = (2/(gamma+1))**(gamma/(gamma-1))
        
        if pressure_ratio > critical_ratio:
            # Unchoked (subsonic) flow
            flow_condition = "Unchoked (Subsonic)"
            term = (gamma/(gamma-1)) * (pressure_ratio**(2/gamma) - pressure_ratio**((gamma+1)/gamma))
            mass_flow_comp = cd_comp * area_comp * np.sqrt(2 * rho1 * p1 * term)
        else:
            # Choked (sonic) flow
            flow_condition = "Choked (Sonic)"
            mass_flow_comp = cd_comp * area_comp * p1 * np.sqrt(gamma / (R_specific * T1) * (2/(gamma+1))**((gamma+1)/(gamma-1)))
        
        st.metric("Flow Condition", flow_condition)
        st.metric("Pressure Ratio (P2/P1)", f"{pressure_ratio:.3f}")
        st.metric("Critical Ratio", f"{critical_ratio:.3f}")
        st.metric("Mass Flow Rate", f"{mass_flow_comp:.4f} kg/s")
        st.metric("Mass Flow Rate", f"{mass_flow_comp*1000:.2f} g/s")
        st.metric("Upstream Density", f"{rho1:.2f} kg/m³")

# Tab 3: Ideal Gas Calculator
with tab3:
    st.header("Ideal Gas Law Calculator")
    st.markdown("Calculate gas properties: PV = mRT or ρ = P/(RT)")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Input Parameters")
        fluid_gas = st.selectbox("Select gas:", 
                                  [f for f in FLUIDS.keys() if FLUIDS[f]['type'] == 'gas'],
                                  key="gas_fluid")
        
        R_spec_gas = FLUIDS[fluid_gas]['R_specific']
        st.info(f"R (specific) = {R_spec_gas} J/(kg·K)")
        
        calc_type = st.radio("Calculate:", ["Density from P,T", "Pressure from ρ,T", "Temperature from P,ρ"])
        
        if calc_type == "Density from P,T":
            p_gas = st.number_input("Pressure (psi):", value=725.0, min_value=0.1, step=10.0, key="p_gas")
            t_gas = st.number_input("Temperature (°C):", value=20.0, step=1.0, key="t_gas")
        elif calc_type == "Pressure from ρ,T":
            rho_gas = st.number_input("Density (kg/m³):", value=50.0, min_value=0.1, step=1.0, key="rho_gas")
            t_gas = st.number_input("Temperature (°C):", value=20.0, step=1.0, key="t_gas2")
        else:  # Temperature from P,ρ
            p_gas = st.number_input("Pressure (psi):", value=725.0, min_value=0.1, step=10.0, key="p_gas2")
            rho_gas = st.number_input("Density (kg/m³):", value=50.0, min_value=0.1, step=1.0, key="rho_gas2")
    
    with col2:
        st.subheader("Results")
        
        if calc_type == "Density from P,T":
            p_pa = p_gas * 6894.76  # psi to Pa
            t_k = t_gas + 273.15
            rho_result = p_pa / (R_spec_gas * t_k)
            st.metric("Density", f"{rho_result:.2f} kg/m³")
            st.metric("Specific Volume", f"{1/rho_result:.4f} m³/kg")
            
        elif calc_type == "Pressure from ρ,T":
            t_k = t_gas + 273.15
            p_result = rho_gas * R_spec_gas * t_k
            st.metric("Pressure", f"{p_result/6894.76:.1f} psi")
            st.metric("Pressure", f"{p_result/1e5:.2f} bar")
            st.metric("Pressure", f"{p_result/1e6:.2f} MPa")
            
        else:  # Temperature from P,ρ
            p_pa = p_gas * 6894.76  # psi to Pa
            t_result = p_pa / (rho_gas * R_spec_gas)
            st.metric("Temperature", f"{t_result:.2f} K")
            st.metric("Temperature", f"{t_result - 273.15:.2f} °C")

# Tab 4: Injector Sizing
with tab4:
    st.header("Injector Orifice Sizing")
    st.markdown("Design orifice diameter for target mass flow rate")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Input Parameters")
        fluid_inj = st.selectbox("Select fluid:", list(FLUIDS.keys()), key="inj_fluid")
        
        target_flow = st.number_input("Target mass flow rate (g/s):", value=100.0, min_value=0.1, step=10.0)
        target_flow_kg = target_flow / 1000  # Convert to kg/s
        
        if FLUIDS[fluid_inj]['type'] == 'liquid':
            use_custom_inj = st.checkbox("Use custom density", key="custom_dens_inj")
            if use_custom_inj:
                density_inj = st.number_input("Density (kg/m³):", value=1000.0, min_value=0.1, key="dens_inj")
            else:
                density_inj = FLUIDS[fluid_inj]['density']
                st.info(f"Using density: {density_inj} kg/m³")
        
        cd_inj = st.number_input("Discharge coefficient (Cd):", value=0.7, min_value=0.1, max_value=1.0, step=0.05, key="cd_inj")
        delta_p_inj = st.number_input("Pressure drop (psi):", value=145.0, min_value=0.1, step=5.0, key="dp_inj")
        
        if FLUIDS[fluid_inj]['type'] == 'liquid':
            num_orifices = st.number_input("Number of orifices:", value=1, min_value=1, step=1)
        
    with col2:
        st.subheader("Results")
        
        if FLUIDS[fluid_inj]['type'] == 'liquid':
            # Incompressible calculation
            # ṁ = Cd * A * ρ * sqrt(2*ΔP/ρ)
            # A = ṁ / (Cd * ρ * sqrt(2*ΔP/ρ))
            
            delta_p_pa_inj = delta_p_inj * 6894.76  # psi to Pa
            flow_per_orifice = target_flow_kg / num_orifices
            
            velocity_inj = np.sqrt(2 * delta_p_pa_inj / density_inj)
            area_required = flow_per_orifice / (cd_inj * density_inj * velocity_inj)
            diameter_required = np.sqrt(4 * area_required / np.pi) / 0.0254  # Convert m to inches
            
            st.metric("Required Diameter (per orifice)", f"{diameter_required:.4f} inches")
            st.metric("Flow per Orifice", f"{flow_per_orifice*1000:.2f} g/s")
            st.metric("Total Area Required", f"{area_required*num_orifices*1550:.4f} in²")
            st.metric("Exit Velocity", f"{velocity_inj:.2f} m/s")
            
        else:
            st.warning("Compressible flow injector sizing requires more parameters. Use the Compressible Flow tab for gas calculations.")

# Tab 5: Mixture Ratio
with tab5:
    st.header("Oxidizer/Fuel Mixture Ratio")
    st.markdown("Calculate mixture ratio and propellant flow rates")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Input Parameters")
        
        calc_mode = st.radio("Calculate:", ["Flow rates from mixture ratio", "Mixture ratio from flow rates"])
        
        oxidizer = st.selectbox("Oxidizer:", [f for f in FLUIDS.keys() if FLUIDS[f]['type'] == 'liquid'], key="ox")
        fuel = st.selectbox("Fuel:", [f for f in FLUIDS.keys() if FLUIDS[f]['type'] == 'liquid'], key="fuel")
        
        if calc_mode == "Flow rates from mixture ratio":
            mixture_ratio = st.number_input("Mixture ratio (O/F):", value=2.5, min_value=0.1, step=0.1)
            total_flow = st.number_input("Total mass flow (g/s):", value=200.0, min_value=0.1, step=10.0)
        else:
            ox_flow = st.number_input("Oxidizer flow (g/s):", value=150.0, min_value=0.0, step=10.0)
            fuel_flow = st.number_input("Fuel flow (g/s):", value=50.0, min_value=0.0, step=10.0)
    
    with col2:
        st.subheader("Results")
        
        if calc_mode == "Flow rates from mixture ratio":
            # ṁ_total = ṁ_ox + ṁ_fuel
            # O/F = ṁ_ox / ṁ_fuel
            # ṁ_ox = (O/F) * ṁ_fuel
            # ṁ_total = (O/F) * ṁ_fuel + ṁ_fuel = ṁ_fuel * (O/F + 1)
            
            fuel_flow_calc = total_flow / (mixture_ratio + 1)
            ox_flow_calc = mixture_ratio * fuel_flow_calc
            
            st.metric("Oxidizer Flow", f"{ox_flow_calc:.2f} g/s")
            st.metric("Fuel Flow", f"{fuel_flow_calc:.2f} g/s")
            st.metric("Mixture Ratio (O/F)", f"{mixture_ratio:.2f}")
            
            # Volume flow rates
            ox_vol = ox_flow_calc / FLUIDS[oxidizer]['density'] * 1000  # L/s
            fuel_vol = fuel_flow_calc / FLUIDS[fuel]['density'] * 1000
            
            st.metric("Oxidizer Volume Flow", f"{ox_vol:.3f} L/s")
            st.metric("Fuel Volume Flow", f"{fuel_vol:.3f} L/s")
            
        else:
            total_flow_calc = ox_flow + fuel_flow
            if fuel_flow > 0:
                mixture_ratio_calc = ox_flow / fuel_flow
                st.metric("Mixture Ratio (O/F)", f"{mixture_ratio_calc:.2f}")
            else:
                st.error("Fuel flow must be > 0")
            
            st.metric("Total Flow", f"{total_flow_calc:.2f} g/s")
            st.metric("Oxidizer Percentage", f"{(ox_flow/total_flow_calc*100):.1f}%")
            st.metric("Fuel Percentage", f"{(fuel_flow/total_flow_calc*100):.1f}%")

st.markdown("---")
st.caption("Built with ❤️ for mission-critical operations | Ishaan Patel")