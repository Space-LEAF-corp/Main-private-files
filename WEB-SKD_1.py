import ui
import {
    toPasskeyTransport,
    toWebAuthnCredential,
} from '@circle-fin/modular-wallets-core'

// 0. retrieve client key and client url from environment vars
const clientKey = import.meta.env.VITE_CLIENT_KEY as string
const clientUrl = import.meta.env.VITE_CLIENT_URL as string

// 1. register or login with a passkey and
//    Create a Passkey Transport from client key
const passkeyTransport = toPasskeyTransport(clientUrl, clientKey)
const credential = await toWebAuthnCredential({
    transport: passkeyTransport,
    mode: WebAuthnMode.Register, //or WebAuthnMode.Login if login
    username: 'your-username'  //replace with actual username
})
v = ui.load_view()
v.present('sheet')