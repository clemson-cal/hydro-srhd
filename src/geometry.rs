use core::ops::Add;
use core::ops::Sub;
use core::ops::Mul;
use core::ops::Div;




// ============================================================================
#[derive(Copy, Clone)]
pub struct Vector3d(pub f64, pub f64, pub f64);

#[derive(Copy, Clone)]
pub enum Direction { X, Y, Z }




// ============================================================================
impl Vector3d
{
    pub fn x(&self) -> f64 { self.0 }
    pub fn y(&self) -> f64 { self.1 }
    pub fn z(&self) -> f64 { self.2 }
}




// ============================================================================
impl Vector3d
{
    pub fn dot(&self, other: &Vector3d) -> f64 {
        self.0 * other.0 + self.1 * other.1 + self.2 * other.2
    }
    pub fn length_squared(&self) -> f64 {
        self.dot(self)
    }
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }
}



// ============================================================================
impl Add<Vector3d> for Vector3d { type Output = Self; fn add(self, u: Vector3d) -> Vector3d { Vector3d(self.0 + u.0, self.1 + u.1, self.2 + u.2) } }
impl Sub<Vector3d> for Vector3d { type Output = Self; fn sub(self, u: Vector3d) -> Vector3d { Vector3d(self.0 - u.0, self.1 - u.1, self.2 - u.2) } }
impl Mul<f64> for Vector3d { type Output = Vector3d; fn mul(self, a: f64) -> Vector3d { Vector3d(self.0 * a, self.1 * a, self.2 * a) } }
impl Div<f64> for Vector3d { type Output = Vector3d; fn div(self, a: f64) -> Vector3d { Vector3d(self.0 / a, self.1 / a, self.2 / a) } }




// ============================================================================
impl Direction
{
    pub fn along(self, other: Direction) -> f64
    {
        match (self, other)
        {
            (Direction::X, Direction::X) => 1.0,
            (Direction::Y, Direction::Y) => 1.0,
            (Direction::Z, Direction::Z) => 1.0,
            _ => 0.0,
        }
    }
    pub fn dot(self, vec: &Vector3d) -> f64
    {
        match self
        {
            Direction::X => vec.x(),
            Direction::Y => vec.y(),
            Direction::Z => vec.z(),
        }           
    }
}




// ============================================================================
impl Mul<f64> for Direction {
    type Output = Vector3d;
    fn mul(self, a: f64) -> Vector3d {
        match self {
            Direction::X => Vector3d(a, 0.0, 0.0),
            Direction::Y => Vector3d(0.0, a, 0.0),
            Direction::Z => Vector3d(0.0, 0.0, a),        
        }
    }
}
