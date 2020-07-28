use num_traits::{One, Zero};
use std::ops::{Add, Mul};
use std::ops::{BitAnd, BitOr, Not};

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Bool<T>(T);
impl From<bool> for Bool<bool> {
    fn from(b: bool) -> Self {
        Bool(b)
    }
}

impl<T> Add for Bool<T>
where
    T: Copy + BitOr<Output = T> + BitAnd + Not,
{
    type Output = Bool<T>;
    fn add(self, o: Self) -> Self::Output {
        Bool(BitOr::bitor(self.0, o.0))
    }
}

impl Zero for Bool<bool> {
    fn zero() -> Self {
        Bool(false)
    }
    fn is_zero(&self) -> bool {
        !self.0
    }
}

impl<T> Mul for Bool<T>
where
    T: Copy + BitOr<Output = T> + BitAnd + Not,
{
    type Output = Self;
    fn mul(self, o: Self) -> Self::Output {
        Bool(BitOr::bitor(self.0, o.0))
    }
}

impl One for Bool<bool> {
    fn one() -> Self {
        Bool(true)
    }
    fn is_one(&self) -> bool {
        self.0
    }
}
