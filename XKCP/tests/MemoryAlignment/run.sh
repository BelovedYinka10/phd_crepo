#!/bin/bash
# Build an XKCP target with clang's or gcc's alignment sanitizer and run the
# misaligned-buffer test against it. Needs compiler-rt installed for clang or
# libubsan installed for gcc.
#
# Usage: run.sh <xkcp-tree> [target]      (default target: generic64)
#
# Pass/fail criteria:
#   - any "runtime error" from UBSan  -> misaligned lane access in the library
#   - nonzero harness exit            -> wrong results on misaligned buffers
set -u
TREE=${1:?usage: run.sh <xkcp-tree> [target]}
TARGET=${2:-generic64}
HERE=$(cd "$(dirname "$0")" && pwd)
COMPILER="${COMPILER:-clang}"
SAN="-fsanitize=alignment -g -O1"

WRAP=$(mktemp /tmp/cc-ubsan.XXXXXX)
printf '#!/bin/sh\nexec %s %s "$@"\n' "$COMPILER" "$SAN" > "$WRAP"
chmod +x "$WRAP"
trap 'rm -f "$WRAP"' EXIT

cd "$TREE"
rm -rf bin/.build
make CC="$WRAP" "$TARGET/libXKCP.so" > /dev/null || { echo "library build failed"; exit 2; }

$COMPILER $SAN -I "bin/$TARGET/libXKCP.so.headers" "$HERE/misalign-test.c" \
    "bin/$TARGET/libXKCP.so" -o "/tmp/misalign-test-$TARGET" || exit 2

LOG=/tmp/misalign-$TARGET.log
LD_LIBRARY_PATH="$TREE/bin/$TARGET" "/tmp/misalign-test-$TARGET" > "$LOG" 2>&1
RC=$?
ERRS=$(grep -c "runtime error" "$LOG")

echo "== $TARGET in $TREE =="
echo "   UBSan alignment errors: $ERRS"
echo "   result mismatches:      $(grep -c MISMATCH "$LOG") (harness rc=$RC)"
grep "runtime error" "$LOG" | sed 's/.*runtime error/   e.g./' | sort -u | head -3
if [ "$ERRS" -eq 0 ] && [ "$RC" -eq 0 ]; then
    echo "   PASS"
    exit 0
else
    echo "   FAIL (full log: $LOG)"
    exit 1
fi
