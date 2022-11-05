set positional-arguments

new version:
    #!/usr/bin/env sh
    if git rev-parse --quiet --verify release/{{version}} > /dev/null; then
        echo "Release branch exists already"
    else
      echo "Creating release branch release/{{version}}"
      git checkout main && \
      git pull && \
      git checkout -b release/{{version}} && \
      sed -i .bck "s/^version =.*$/version = \"{{version}}\"/" ./Cargo.toml && \
      cargo check && \
      git commit -am "Prepare release branch {{version}}" && \
      git push -u origin release/{{version}}
    fi

@check version:
    git checkout release/{{version}} && git pull
    echo "Running linter and unittests"
    cargo clippy && cargo fmt && cargo test -q && cargo doc

@release version:
    git tag {{version}}
    git push --tags
    cargo publish

test:
    #!/usr/bin/env zsh
    echo -ne "Checking formatting and doc generation"
    (cargo clippy && cargo fmt --check && cargo test -q && cargo doc && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

coverage:
    #!/usr/bin/env zsh
    export RUSTFLAGS="-Cinstrument-coverage"
    export CARGO_INCREMENTAL=0
    export RUSTFLAGS="-Zprofile -Ccodegen-units=1 -Copt-level=0 -Clink-dead-code -Coverflow-checks=off -Zpanic_abort_tests -Cpanic=abort"
    export RUSTDOCFLAGS="-Cpanic=abort"
    cargo +nightly build
    export LLVM_PROFILE_FILE="atg-%p-%m.profraw"
    cargo +nightly test
    grcov . -s . --binary-path ./target/debug/ -t html --branch --ignore-not-existing -o ./target/debug/coverage/
    open target/debug/coverage/index.html
