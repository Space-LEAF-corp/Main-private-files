#![no_std]
#![no_main]

mod turbo_stack;
mod shadow_panther;
mod raptor_core;
mod config;

use core::panic::PanicInfo;

#[no_mangle]
pub extern "C" fn _start() -> ! {
    // Boot sequence for Raptor OS
    config::init();
    turbo_stack::init_turbo_stack();
    shadow_panther::init_shadow_panther_copilot();
    raptor_core::run_shell();

    loop {}
}

#[panic_handler]
fn panic(_info: &PanicInfo) -> ! {
    loop {}
}
