ethics_config = EthicsConfig()
ethics_engine = EthicsEngine(ethics_config)

def handle_client(client):
    username = authenticate(client)
    if not username:
        return

    notify_user(client, "Welcome to the private ethical chat. Use 'bullshit' wisely to flag issues.")

    kicked = False

    while True:
        try:
            message = client.recv(1024).decode('utf-8').strip()
            if not message:
                break

            timestamp = datetime.now()
            CHAT_HISTORY.append((username, message, timestamp))
            print(f"[{timestamp}] {username}: {message}")

            actions = ethics_engine.evaluate_message(username, message, CHAT_HISTORY)

            # Apply actions
            for act in actions:
                if act.action == ActionType.NONE:
                    notify_user(client, act.reason)
                elif act.action == ActionType.WARN:
                    notify_user(client, act.reason)
                elif act.action == ActionType.COOLDOWN:
                    notify_user(client, act.reason)
                elif act.action == ActionType.FLAG:
                    notify_user(client, act.reason)
                    # TODO: notify moderators / log
                elif act.action == ActionType.BAN:
                    notify_user(client, act.reason)
                    kicked = True
                    client.close()
                    if username in AUTHENTICATED_USERS:
                        AUTHENTICATED_USERS.remove(username)
                    break

            if kicked:
                break

        except:
            break

    if not kicked:
        kick_user(client, username, "Disconnected.")