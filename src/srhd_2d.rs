use std::ops::{Add, Sub, Mul, Div};
use std::cell::Cell;
use crate::geometry::{Direction, Vector3d};




thread_local! {
    static PRESSURE_GUESS: Cell<Option<f64>> = Cell::new(Some(0.0));
}




// ============================================================================
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature="hdf5", repr(C), derive(hdf5::H5Type))]
#[cfg_attr(feature="serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Conserved(pub f64, pub f64, pub f64, pub f64);

#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature="hdf5", repr(C), derive(hdf5::H5Type))]
#[cfg_attr(feature="serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Primitive(pub f64, pub f64, pub f64, pub f64);




// ============================================================================
impl Add<Primitive> for Primitive { type Output = Self; fn add(self, u: Primitive) -> Primitive { Primitive(self.0 + u.0, self.1 + u.1, self.2 + u.2, self.3 + u.3) } }
impl Sub<Primitive> for Primitive { type Output = Self; fn sub(self, u: Primitive) -> Primitive { Primitive(self.0 - u.0, self.1 - u.1, self.2 - u.2, self.3 - u.3) } }
impl Mul<f64> for Primitive { type Output = Primitive; fn mul(self, a: f64) -> Primitive { Primitive(self.0 * a, self.1 * a, self.2 * a, self.3 * a) } }
impl Div<f64> for Primitive { type Output = Primitive; fn div(self, a: f64) -> Primitive { Primitive(self.0 / a, self.1 / a, self.2 / a, self.3 / a) } }




// ============================================================================
impl Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1, self.2 + u.2, self.3 + u.3) } }
impl Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1, self.2 - u.2, self.3 - u.3) } }
impl Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a, self.2 * a, self.3 * a) } }
impl Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a, self.2 / a, self.3 / a) } }




// ============================================================================
impl From<Primitive> for [f64; 4] { fn from(a: Primitive) -> [f64; 4] { [a.0, a.1, a.2, a.3] } }
impl From<Conserved> for [f64; 4] { fn from(a: Conserved) -> [f64; 4] { [a.0, a.1, a.2, a.3] } }
impl From<[f64; 4]> for Primitive { fn from(a: [f64; 4]) -> Primitive { Primitive(a[0], a[1], a[2], a[3]) } }
impl From<[f64; 4]> for Conserved { fn from(a: [f64; 4]) -> Conserved { Conserved(a[0], a[1], a[2], a[3]) } }

impl Default for Conserved { fn default() -> Self { Conserved(0.0, 0.0, 0.0, 0.0) } }
impl Default for Primitive { fn default() -> Self { Primitive(0.0, 0.0, 0.0, 0.0) } }

impl Conserved { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e } }
impl Primitive { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e } }

pub enum RecoveredPrimitive {
    Success(Primitive),
    NegativePressure(Primitive),
    RootFinderFailed(Conserved),
}

impl RecoveredPrimitive {
    pub fn unwrap(self) -> Primitive {
        match self {
            Self::Success(p) => p,
            Self::NegativePressure(p) => panic!("primitive recovery failed with negative pressure {:?}", p),
            Self::RootFinderFailed(u) => panic!("primitive recovery failed on conserved state {:?}", u),
        }
    }
}




// ============================================================================
impl Conserved {
    pub fn lab_frame_density(self)  -> f64 { self.0 }
    pub fn momentum_1       (self)  -> f64 { self.1 }
    pub fn momentum_2       (self)  -> f64 { self.2 }
    pub fn energy_density   (self)  -> f64 { self.3 }

    pub fn momentum_vector(self) -> Vector3d {
        Vector3d(self.momentum_1(), self.momentum_2(), 0.0)
    }

    pub fn momentum(self, direction: Direction) -> f64 {
        match direction {
            Direction::X => self.momentum_1(),
            Direction::Y => self.momentum_2(),
            Direction::Z => 0.0,
        }
    }

    pub fn momentum_squared(self) -> f64 {
        let s1 = self.momentum_1();
        let s2 = self.momentum_2();
        return s1 * s1 + s2 * s2;
    }

    pub fn to_primitive(self, gamma_law_index: f64) -> RecoveredPrimitive {
        use RecoveredPrimitive::*;

        PRESSURE_GUESS.with(|pressure_guess_cell| {

            if let Some(pressure_guess) = pressure_guess_cell.get() {

                // The pressure guess is Some, so we are in normal mode. If
                // the inversion succeeds, set the guess pressure to the new
                // result pressure, and return the primitive state. If the
                // inversion fails, then set the pressure guess to None to
                // indicate safety mode, and recurse.

                match self.to_primitive_internal(gamma_law_index, pressure_guess) {
                    Success(prim) => {
                        pressure_guess_cell.set(Some(prim.gas_pressure()));
                        Success(prim)
                    }
                    _ => {
                        pressure_guess_cell.set(None);
                        self.to_primitive(gamma_law_index)
                    }
                }
            } else {

                // The pressure guess is None, so we are in safety mode. If
                // the inversion succeeds, set the guess pressure to the new
                // result pressure, and return the primitive state. If the
                // inversion fails, then restore the state to normal mode and
                // return the error state.

                match self.to_primitive_internal(gamma_law_index, 0.0) {
                    Success(prim) => {
                        pressure_guess_cell.set(Some(prim.gas_pressure()));
                        Success(prim)
                    }
                    NegativePressure(p) => {
                        pressure_guess_cell.set(Some(0.0));
                        NegativePressure(p)
                    }
                    RootFinderFailed(u) => {
                        pressure_guess_cell.set(Some(0.0));
                        RootFinderFailed(u)
                    }
                }
            }
        })
    }

    fn to_primitive_internal(self, gamma_law_index: f64, pressure_guess: f64) -> RecoveredPrimitive {

        let newton_iter_max = 50;
        let error_tolerance = 1e-12 * self.lab_frame_density();
        let gm              = gamma_law_index;
        let m               = self.lab_frame_density();
        let tau             = self.energy_density();
        let ss              = self.momentum_squared();
        let w0;
        let mut iteration   = 0;
        let mut p           = pressure_guess;

        loop {
            let et = tau + p + m;
            let b2 = f64::min(ss / et / et, 1.0 - 1e-10);
            let w2 = 1.0 / (1.0 - b2);
            let w  = f64::sqrt(w2);
            let e  = (tau + m * (1.0 - w) + p * (1.0 - w2)) / (m * w);
            let d  = m / w;
            let h  = 1.0 + e + p / d;
            let a2 = gm * p / (d * h);
            let f  = d * e * (gm - 1.0) - p;
            let g  = b2 * a2 - 1.0;

            p -= f / g;

            if f64::abs(f) < error_tolerance || iteration == newton_iter_max {
                w0 = w;
                break;
            }
            iteration += 1;
        }

        let result = Primitive(
            m / w0,
            w0 * self.momentum_1() / (tau + m + p),
            w0 * self.momentum_2() / (tau + m + p),
            p
        );

        if iteration == newton_iter_max {
            RecoveredPrimitive::RootFinderFailed(self)
        } else if p < 0.0 {
            RecoveredPrimitive::NegativePressure(result)
        } else {
            RecoveredPrimitive::Success(result)
        }
    }
}




// ============================================================================
impl Primitive {
    pub fn mass_density(self) -> f64 { self.0 }
    pub fn gamma_beta_1(self) -> f64 { self.1 }
    pub fn gamma_beta_2(self) -> f64 { self.2 }
    pub fn gas_pressure(self) -> f64 { self.3 }

    pub fn velocity_1(self) -> f64 {
        return self.gamma_beta_1() / self.lorentz_factor();
    }

    pub fn velocity_2(self) -> f64 {
        return self.gamma_beta_2() / self.lorentz_factor();
    }

    pub fn gamma_beta(self, direction: Direction) -> f64 {
        match direction {
            Direction::X => self.gamma_beta_1(),
            Direction::Y => self.gamma_beta_2(),
            Direction::Z => 0.0,
        }
    }

    pub fn velocity(self, direction: Direction) -> f64 {
        match direction {
            Direction::X => self.velocity_1(),
            Direction::Y => self.velocity_2(),
            Direction::Z => 0.0,
        }
    }

    pub fn four_velocity_vector(self) -> Vector3d {
        Vector3d(self.gamma_beta_1(), self.gamma_beta_2(), 0.0) / self.lorentz_factor()
    }

    pub fn gamma_beta_squared(self) -> f64 {
        let u1 = self.gamma_beta_1();
        let u2 = self.gamma_beta_2();
        return u1 * u1 + u2 * u2;
    }

    pub fn lorentz_factor_squared(self) -> f64 {
        return 1.0 + self.gamma_beta_squared();
    }

    pub fn lorentz_factor(self) -> f64 {
        return f64::sqrt(self.lorentz_factor_squared());
    }

    pub fn enthalpy_density(self, gamma_law_index: f64) -> f64 {
        return self.mass_density() + self.gas_pressure() * (1.0 + 1.0 / (gamma_law_index - 1.0));
    }

    pub fn specific_enthalpy(self, gamma_law_index: f64) -> f64 {
        return self.enthalpy_density(gamma_law_index) / self.mass_density();
    }

    pub fn specific_entropy(self, gamma_law_index: f64) -> f64 {
        return f64::ln(self.gas_pressure() / f64::powf(self.mass_density(), gamma_law_index));
    }

    pub fn sound_speed_squared(self, gamma_law_index: f64) -> f64 {
        return gamma_law_index * self.gas_pressure() / self.enthalpy_density(gamma_law_index);
    }

    pub fn outer_wavespeeds(self, direction: Direction, gamma_law_index: f64) -> (f64, f64) {
        let a2 = self.sound_speed_squared(gamma_law_index);
        let uu = self.gamma_beta_squared();
        let vn = self.velocity(direction);
        let vv = uu / (1.0 + uu);
        let v2 = vn * vn;
        let k0 = f64::sqrt(a2 * (1.0 - vv) * (1.0 - vv * a2 - v2 * (1.0 - a2)));

        return ((vn * (1.0 - a2) - k0) / (1.0 - vv * a2),
                (vn * (1.0 - a2) + k0) / (1.0 - vv * a2));
    }

    pub fn max_signal_speed(self, gamma_law_index: f64) -> f64 {
        f64::sqrt(self.gamma_beta_squared() + self.sound_speed_squared(gamma_law_index))
    }

    pub fn to_conserved(self, gamma_law_index: f64) -> Conserved {
        let w = self.lorentz_factor();
        let m = self.mass_density() * w;
        let h = self.specific_enthalpy(gamma_law_index);

        return Conserved(
            m,
            m * h * self.gamma_beta_1(),
            m * h * self.gamma_beta_2(),
            m * (h * w - 1.0) - self.gas_pressure());
    }

    pub fn flux_vector(self, direction: Direction, gamma_law_index: f64) -> Conserved {
        let pg = self.gas_pressure();
        let vn = self.velocity(direction);

        let pressure_term = Conserved(
            0.0,
            pg * direction.along(Direction::X),
            pg * direction.along(Direction::Y),
            pg * vn);

        let advective_term = self.to_conserved(gamma_law_index) * vn;

        return advective_term + pressure_term;
    }
}




// ============================================================================
impl Primitive
{
    pub fn spherical_geometry_source_terms(self, spherical_radius: f64, polar_angle_theta: f64, gamma_law_index: f64) -> Conserved
    {
        let cotq = f64::tan(std::f64::consts::FRAC_PI_2 - polar_angle_theta);
        let ur = self.gamma_beta_1();
        let uq = self.gamma_beta_2();
        let up = 0.0;
        let pg = self.gas_pressure();
        let h0 = self.enthalpy_density(gamma_law_index);
        let sd = 0.0;
        let sr = (2.0  * pg + h0 * (uq * uq        + up * up)) / spherical_radius;
        let sq = (cotq * pg + h0 * (up * up * cotq - ur * uq)) / spherical_radius;
        // let sp =        -up * h0 * (ur + uq * cotq) / spherical_radius;
        let se = 0.0;
        Conserved(sd, sr, sq, se)
    }
}




// ============================================================================
pub fn riemann_hlle(pl: Primitive, pr: Primitive, direction: Direction, gamma_law_index: f64) -> Conserved {
    let ul = pl.to_conserved(gamma_law_index);
    let ur = pr.to_conserved(gamma_law_index);
    let fl = pl.flux_vector(direction, gamma_law_index);
    let fr = pr.flux_vector(direction, gamma_law_index);

    let (alm, alp) = pl.outer_wavespeeds(direction, gamma_law_index);
    let (arm, arp) = pr.outer_wavespeeds(direction, gamma_law_index);
    let ap = alp.max(arp).max(0.0);
    let am = alm.min(arm).min(0.0);

    (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am)
}




// ============================================================================
pub enum RiemannSolverMode
{
    HlleFlux,
    HlleFluxAcrossMovingFace(f64),
    HllcFlux,
    HllcFluxAcrossMovingFace(f64),
    HllcFluxAcrossContact,
}




// ============================================================================
pub fn riemann_hllc(pl: Primitive, pr: Primitive, nhat: Direction, gamma_law_index: f64, mode: RiemannSolverMode) -> (Conserved, f64)
{
    let ul = pl.to_conserved(gamma_law_index);
    let ur = pr.to_conserved(gamma_law_index);
    let fl = pl.flux_vector(nhat, gamma_law_index);
    let fr = pr.flux_vector(nhat, gamma_law_index);

    let (alm, alp) = pl.outer_wavespeeds(nhat, gamma_law_index);
    let (arm, arp) = pr.outer_wavespeeds(nhat, gamma_law_index);
    let ar = f64::max(alp, arp);
    let al = f64::min(alm, arm);

    // Equations (9) and (11)
    let u_hll = (ur * ar - ul * al + (fl - fr))           / (ar - al);
    let f_hll = (fl * ar - fr * al - (ul - ur) * ar * al) / (ar - al);

    let discriminant = |a: f64, b: f64, c: f64| -> f64
    {
        b * b - 4.0 * a * c
    };

    let quadratic_root = |a: f64, b: f64, c: f64| -> f64
    {
        let d = discriminant(a, b, c);
        if d < 0.0 {
            0.0
        } else if a.abs() < 1e-8 {
            -c / b
        } else {
            (-b - d.sqrt()) / 2.0 / a
        }
    };

    // Equation (18) for a-star and p-star
    let a_star_and_p_star = || -> (f64, f64)
    {
        // Mignone defines total energy to include rest mass
        let ue_hll = u_hll.energy_density() + u_hll.lab_frame_density();
        let fe_hll = f_hll.energy_density() + f_hll.lab_frame_density();
        let um_hll = u_hll.momentum(nhat);
        let fm_hll = f_hll.momentum(nhat);
        let a_star = quadratic_root(fe_hll, -fm_hll - ue_hll, um_hll);
        let p_star = -fe_hll * a_star + fm_hll;
        if a_star.is_nan() {
            panic!("a* is NaN {:?} {:?}", pl, pr);
        }
        (a_star, p_star)
    };

    // Equations (16)
    let star_state_flux = |u: Conserved, f: Conserved, p: Primitive, a: f64, vface: f64, a_star: f64, p_star: f64| -> Conserved
    {
        let e = u.energy_density() + u.lab_frame_density();
        let s = u.momentum_vector();
        let m = u.momentum(nhat);
        let v = p.velocity(nhat);
        let n = s - nhat * m; // transverse momentum
        let es = (e * (a - v) + p_star * a_star - p.gas_pressure() * v) / (a - a_star);
        let ms = (m * (a - v) + p_star          - p.gas_pressure())     / (a - a_star);
        let ns =  n * (a - v)                                           / (a - a_star);
        let ds = u.lab_frame_density() * (a - v)                        / (a - a_star);
        let ss = nhat * ms + ns;
        let us = Conserved(ds, ss.x(), ss.y(), es - ds);
        let fs = f + (us - u) * a;
        fs - us * vface
    };

    match mode
    {
        RiemannSolverMode::HlleFlux => {
            if      0.0 < al { (fl   , 0.0) }
            else if 0.0 > ar { (fr   , 0.0) }
            else             { (f_hll, 0.0) }
        }

        RiemannSolverMode::HlleFluxAcrossMovingFace(vface) => {
            if      vface < al { (fl    - ul    * vface, vface) }
            else if vface > ar { (fr    - ur    * vface, vface) }
            else               { (f_hll - u_hll * vface, vface) }
        }

        RiemannSolverMode::HllcFlux => {
            let (a_star, p_star) = a_star_and_p_star();
            if      0.0 <= al     { (fl, 0.0) }
            else if 0.0 >= ar     { (fr, 0.0) }
            else if 0.0 <= a_star { (star_state_flux(ul, fl, pl, al, 0.0, a_star, p_star), 0.0) }
            else if 0.0 >= a_star { (star_state_flux(ur, fr, pr, ar, 0.0, a_star, p_star), 0.0) }
            else                  { unreachable!() }
        }

        RiemannSolverMode::HllcFluxAcrossMovingFace(vface) => {
            let (a_star, p_star) = a_star_and_p_star();
            if      vface <= al     { (fl - ul * vface, vface) }
            else if vface >= ar     { (fr - ur * vface, vface) }
            else if vface <= a_star { (star_state_flux(ul, fl, pl, al, vface, a_star, p_star), vface) }
            else if vface >= a_star { (star_state_flux(ur, fr, pr, ar, vface, a_star, p_star), vface) }
            else                    { unreachable!() }
        }

        RiemannSolverMode::HllcFluxAcrossContact => {
            let (a_star, p_star) = a_star_and_p_star();
            (star_state_flux(ul, fl, pl, al, a_star, a_star, p_star), a_star)
        }
    }
}




/**
 * Return HLL[E/C] Godunov fluxes of the conserved quantities and a passive
 * scalar. The inputs sl and sr are primitive-like: they are the comoving scalar
 * density, just like the primitives contain the comoving mass density.
 */
pub fn riemann_hllc_scalar(
    pl: Primitive,
    pr: Primitive,
    sl: f64,
    sr: f64,
    nhat: Direction,
    gamma_law_index: f64,
    mode: RiemannSolverMode) -> (Conserved, f64, f64)
{
    let cl = sl / pl.mass_density(); // scalar concentration
    let cr = sr / pr.mass_density();
    let (f, v) = riemann_hllc(pl, pr, nhat, gamma_law_index, mode);

    let g = if f.lab_frame_density() < 0.0 {
        cr * f.lab_frame_density()
    } else {
        cl * f.lab_frame_density()
    };
    (f, g, v)
}




// ============================================================================
#[cfg(test)]
mod tests
{
    use super::*;

    fn panic_unless_recovery_is_accurate(primitive: Primitive)
    {
        let gamma_law_index = 4.0 / 3.0;
        let u = primitive.to_conserved(gamma_law_index);
        let p = u.to_primitive(gamma_law_index).unwrap();
        println!("{:?} {:?}", primitive, p);
        assert!(f64::abs(primitive.gamma_beta_squared() - p.gamma_beta_squared()) < 1e-10);
        assert!(f64::abs(primitive.gas_pressure() - p.gas_pressure()) / primitive.gas_pressure() < 1e-10);
    }

    #[test]
    fn can_recover_primitive()
    {
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.2, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 0.2, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.5, 0.5, 1e-3));
        panic_unless_recovery_is_accurate(Primitive(1.0, 5.0, 5.0, 1e+3));
    }

    #[test]
    fn can_recover_primitive_hard()
    {
        panic_unless_recovery_is_accurate(Primitive(1e+12, 0.0, 0.0, 1e+12));        
        panic_unless_recovery_is_accurate(Primitive(1e-12, 0.0, 0.0, 1e-12));
        panic_unless_recovery_is_accurate(Primitive(1e+12, 0.5, 0.5, 1e+12));
        panic_unless_recovery_is_accurate(Primitive(1e-12, 0.5, 0.5, 1e-12));
        panic_unless_recovery_is_accurate(Primitive(1e+16, 0.5, 0.5, 1e+16));
        panic_unless_recovery_is_accurate(Primitive(1e-16, 0.5, 0.5, 1e-16));
        panic_unless_recovery_is_accurate(Primitive(1e+26, 0.5, 0.5, 1e+26));
        panic_unless_recovery_is_accurate(Primitive(1e-26, 0.5, 0.5, 1e-26));
    }

    #[test]
    fn can_obtain_hlle_flux()
    {
        let p = Primitive(1.0, 0.0, 0.0, 1.0);
        let f = riemann_hlle(p, p, Direction::X, 4.0 / 3.0);
        assert!((f - p.flux_vector(Direction::X, 4.0 / 3.0)).small(1e-12));
    }

    #[test]
    fn can_obtain_hllc_flux_x()
    {
        let pl = Primitive(1.0, 0.5, 0.0, 1.0);
        let pr = Primitive(0.1, 0.5, 0.0, 1.0);
        let (fstar, vface) = riemann_hllc(pl, pr, Direction::X, 4.0 / 3.0, RiemannSolverMode::HllcFluxAcrossContact);

        assert_eq!(fstar.lab_frame_density(), 0.0);
        assert_eq!(vface, pl.velocity_1());
    }

    #[test]
    fn can_obtain_hllc_flux_y()
    {
        let pl = Primitive(1.0, 0.0, 0.5, 1.0);
        let pr = Primitive(0.1, 0.0, 0.5, 1.0);
        let (fstar, vface) = riemann_hllc(pl, pr, Direction::Y, 4.0 / 3.0, RiemannSolverMode::HllcFluxAcrossContact);

        assert_eq!(fstar.lab_frame_density(), 0.0);
        assert_eq!(vface, pl.velocity_2());
    }
}
