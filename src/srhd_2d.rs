#![allow(unused)]

// ============================================================================
#[derive(Copy, Clone)]
pub struct Conserved(pub f64, pub f64, pub f64, pub f64);

#[derive(Copy, Clone)]
pub struct Primitive(pub f64, pub f64, pub f64, pub f64);

#[derive(Copy, Clone)]
pub enum Direction { X, Y }




// ============================================================================
impl Direction
{
    fn dot(self, other: Direction) -> f64
    {
        match (self, other)
        {
            (Direction::X, Direction::X) => 1.0,
            (Direction::Y, Direction::Y) => 1.0,
            _ => 0.0,
        }
    }
}




// ============================================================================
impl std::ops::Add<Primitive> for Primitive { type Output = Self; fn add(self, u: Primitive) -> Primitive { Primitive(self.0 + u.0, self.1 + u.1, self.2 + u.2, self.3 + u.3) } }
impl std::ops::Sub<Primitive> for Primitive { type Output = Self; fn sub(self, u: Primitive) -> Primitive { Primitive(self.0 - u.0, self.1 - u.1, self.2 - u.2, self.3 - u.3) } }
impl std::ops::Mul<f64> for Primitive { type Output = Primitive; fn mul(self, a: f64) -> Primitive { Primitive(self.0 * a, self.1 * a, self.2 * a, self.3 * a) } }
impl std::ops::Div<f64> for Primitive { type Output = Primitive; fn div(self, a: f64) -> Primitive { Primitive(self.0 / a, self.1 / a, self.2 / a, self.3 / a) } }




// ============================================================================
impl std::ops::Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1, self.2 + u.2, self.3 + u.3) } }
impl std::ops::Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1, self.2 - u.2, self.3 - u.3) } }
impl std::ops::Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a, self.2 * a, self.3 * a) } }
impl std::ops::Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a, self.2 / a, self.3 / a) } }




// ============================================================================
impl Into<[f64; 4]> for Primitive {
    fn into(self) -> [f64; 4] {
        [self.0, self.1, self.2, self.3]
    }    
}

impl From<[f64; 4]> for Primitive {
    fn from(a:  [f64; 4]) -> Primitive {
        Primitive(a[0], a[1], a[2], a[3])
    }
}

impl Conserved { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e } }
impl Primitive { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e } }




// ============================================================================
impl Conserved {
    pub fn lab_frame_density(self)  -> f64 { self.0 }
    pub fn momentum_1       (self)  -> f64 { self.1 }
    pub fn momentum_2       (self)  -> f64 { self.2 }
    pub fn energy_density   (self)  -> f64 { self.3 }

    pub fn momentum_squared(self) -> f64 {
        let s1 = self.momentum_1();
        let s2 = self.momentum_2();
        return s1 * s1 + s2 * s2;
    }

    pub fn to_primitive(self, gamma_law_index: f64) -> Primitive {
        let newton_iter_max = 50;
        let error_tolerance = 1e-12;
        let gm              = gamma_law_index;
        let m               = self.lab_frame_density();
        let tau             = self.energy_density();
        let ss              = self.momentum_squared();
        let mut found       = false;
        let mut iteration   = 0;
        let mut w0          = 1.0;
        let mut p           = 0.0;

        while iteration <= newton_iter_max {
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

            if f64::abs(f) < error_tolerance || (f64::abs(f) < error_tolerance && iteration == newton_iter_max) {
                w0 = w;
                found = true;
                break;
            }
            iteration += 1;
        }

        return Primitive(
            m / w0,
            w0 * self.momentum_1() / (tau + m + p),
            w0 * self.momentum_2() / (tau + m + p),
            p
        );
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
        }
    }

    pub fn velocity(self, direction: Direction) -> f64 {
        match direction {
            Direction::X => self.velocity_1(),
            Direction::Y => self.velocity_2(),
        }
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
            pg * direction.dot(Direction::X),
            pg * direction.dot(Direction::Y),
            pg * vn);

        let advective_term = self.to_conserved(gamma_law_index) * vn;

        return advective_term + pressure_term;
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
#[cfg(test)]
mod tests
{
    use super::*;

    fn panic_unless_recovery_is_accurate(primitive: Primitive)
    {
        let gamma_law_index = 4.0 / 3.0;
        let u = primitive.to_conserved(gamma_law_index);
        let p = u.to_primitive(gamma_law_index);
        assert!((primitive - p).small(1e-10));
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
}
