# polymurhash-zig

Zig port of [PolymurHash](https://github.com/orlp/polymur-hash) by [Orson Peters](https://github.com/orlp).

## Usage

1. Create or modify the `build.zig.zon` file in the project root to include `polymurhash-zig` as a dependency.
    
    <details>

    <summary><code>build.zig.zon</code> example</summary>

    ```zig
    .{
        .name = "<name of your program>",
        .version = "<version of your program>",
        .dependencies = .{
            .polymurhash = .{
                .url = "https://github.com/e4m2/polymurhash-zig/archive/refs/tags/<git tag>.tar.gz",
                .hash = "<package hash>",
            },
        },
    }
    ```

    If unsure what to fill out for `<package hash>`, set it to `12200000000000000000000000000000000000000000000000000000000000000000` and Zig will tell you the correct value in an error message.

    </details>

2. Add `polymurhash-zig` as a dependency in `build.zig`.

    <details>

    <summary><code>build.zig</code> example</summary>

    ```zig
    const polymurhash = b.dependency("polymurhash", .{});
    exe.addModule("polymurhash", polymurhash.module("polymurhash"));
    ```

    </details>