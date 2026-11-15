# ETH VPN from the macOS terminal, one command, no typing

This sets up one command, `ethvpn`, that connects to the ETH VPN without the Cisco
app. Your password is stored in the macOS Keychain and the 6-digit 2FA code is
generated on your Mac, so connecting asks for nothing.

Throughout, replace:
- `<user>` with your ETH username;
- `<group>` with your VPN group: `student-net`, `staff-net`, or a department
  group (for example `bsse-internal`). Ask your IT support which one gives
  access to your servers.

## 1. Install

Needs [Homebrew](https://brew.sh) and the default zsh shell.

```sh
brew install openconnect oath-toolkit
```

If the Cisco Secure Client app is running, disconnect and quit it. Only one VPN
client can be connected at a time.

## 2. Get your 2FA secret

Your authenticator app generates codes from a secret key: about 32 characters
of A–Z and 2–7. You need that key.

- **Easiest:** reset your ETH MFA at password.ethz.ch and re-register. When the
  QR code appears, look for an option to show the key as text (often "can't
  scan?"), copy it, then scan the QR with your phone app as usual.
- **From Google Authenticator** (Transfer accounts → Export): save a screenshot
  of the QR code as `qr.png`. Then follow the steps in the appendix. Export
  only the ETH account, so it fits in one QR code.

## 3. Save password and secret in the Keychain

```sh
read -s "pw?ETH network password: "; echo
security add-generic-password -U -a <user> -s eth_radius -w "$pw"; unset pw
read -s "sec?2FA secret key: "; echo
security add-generic-password -U -a <user> -s eth_otp_secret -w "$sec"; unset sec
```

Check it: this should print the same code your phone shows.

```sh
oathtool --totp=sha1 -b "$(security find-generic-password -a <user> -s eth_otp_secret -w)"
```

## 4. Add the command to `~/.zshrc`

```sh
ethvpn() {
  local pw sec
  pw=$(security find-generic-password -a <user> -s eth_radius -w)
  sec=$(security find-generic-password -a <user> -s eth_otp_secret -w)
  [ -n "$pw" ] && [ -n "$sec" ] || { echo "Keychain entries missing" >&2; return 1; }
  printf '%s\n%s\n' "$pw" "$(oathtool --totp=sha1 -b "$sec")" |
    sudo openconnect --passwd-on-stdin --user='<user>@<group>.ethz.ch' --usergroup=<group> \
      --useragent=AnyConnect --no-external-auth sslvpn.ethz.ch
}
```

Then run `source ~/.zshrc`.

**Optional: skip the sudo password.** Run `sudo visudo -f /etc/sudoers.d/openconnect`
and add this line. Use your Mac username, and the path that `which openconnect`
prints (`/usr/local/bin/openconnect` on Intel Macs):

```
<mac-user> ALL=(root) NOPASSWD: /opt/homebrew/bin/openconnect
```

## Use

```sh
ethvpn      # connects; Ctrl-C disconnects
```

- You can ignore a final `** Error: The parameters were not valid` line; the VPN works.
- If you get `Login failed`, check three things: the password (redo step 3),
  whether the step 3 check matches your phone (if not, get the secret again),
  and whether `<group>` is a group your account may use.
- If your ETH password changes, redo step 3.

## Appendix: secret from a Google Authenticator export

The export QR contains an `otpauth-migration://` link. Read it with the macOS QR
reader (`zbarimg` often fails on these). The `swift` command needs the Xcode command-line
tools; if it's missing, run `xcode-select --install`.

```sh
cat > qr.swift <<'EOF'
import CoreImage
let img = CIImage(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1]))!
let det = CIDetector(ofType: CIDetectorTypeQRCode, context: nil, options: [CIDetectorAccuracy: CIDetectorAccuracyHigh])!
for case let f as CIQRCodeFeature in det.features(in: img) { print(f.messageString ?? "") }
EOF
swift qr.swift qr.png > qr.txt
```

Extract the secrets. Each line shows an index, the account name, the issuer and the secret:

```sh
python3 - qr.txt <<'EOF'
import sys, base64, urllib.parse
def fields(b, i=0):
    while i < len(b):
        k, i = varint(b, i)
        if k & 7 == 0: v, i = varint(b, i)
        else: n, i = varint(b, i); v, i = b[i:i+n], i+n
        yield k >> 3, v
def varint(b, i, r=0, s=0):
    while True:
        c = b[i]; i += 1; r |= (c & 0x7f) << s; s += 7
        if c < 0x80: return r, i
q = urllib.parse.urlparse(open(sys.argv[1]).read().strip()).query
accts = [v for f, v in fields(base64.b64decode(urllib.parse.parse_qs(q)['data'][0])) if f == 1]
for n, v in enumerate(accts):
    d = dict(fields(v))
    print(n, d.get(2, b'').decode(), d.get(3, b'').decode(), base64.b32encode(d[1]).decode().rstrip('='))
EOF
```

Copy the secret from the ETH line into step 3. Then delete the files:
`rm qr.png qr.txt qr.swift`. Anyone who has them can generate your codes.
