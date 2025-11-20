.PHONY: build release test readme

build:
	cargo build
	cargo run --bin update_readme

release:
	cargo build --release
	cargo run --bin update_readme

test:
	cargo test

readme:
	cargo run --bin update_readme