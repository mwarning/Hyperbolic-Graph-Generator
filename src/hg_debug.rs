
// The debug version
#[cfg(feature = "debug")]
#[macro_export]
macro_rules! hg_debug {
    ($( $args:expr ),*) => { eprintln!( $( $args ),* ); }
}

// Non-debug version
#[cfg(not(feature = "debug"))]
#[macro_export]
macro_rules! hg_debug {
    ($( $args:expr ),*) => {}
}
