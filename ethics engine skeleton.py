class EthicsEngine:
    def __init__(self, config):
        self.config = config
        self.banned_users = {}      # username -> datetime
        self.cooldowns = {}         # username -> datetime

    def is_banned(self, username: str) -> bool:
        now = datetime.now()
        if username in self.banned_users and self.banned_users[username] > now:
            return True
        if username in self.banned_users and self.banned_users[username] <= now:
            del self.banned_users[username]
        return False

    def is_on_cooldown(self, username: str) -> bool:
        now = datetime.now()
        if username in self.cooldowns and self.cooldowns[username] > now:
            return True
        if username in self.cooldowns and self.cooldowns[username] <= now:
            del self.cooldowns[username]
        return False

    def record_ban(self, username: str, duration: timedelta) -> ModerationAction:
        until = datetime.now() + duration
        self.banned_users[username] = until
        return ModerationAction(
            action=ActionType.BAN,
            reason=f"Banned until {until}",
            ban_until=until
        )

    def record_cooldown(self, username: str, duration: timedelta) -> ModerationAction:
        until = datetime.now() + duration
        self.cooldowns[username] = until
        return ModerationAction(
            action=ActionType.COOLDOWN,
            reason=f"Cooldown until {until}",
            cooldown_until=until
        )

    def evaluate_message(self, username: str, message: str, chat_history: list) -> list[ModerationAction]:
        """
        Returns a list of ModerationAction objects.
        This is where your JARVONDIS ethic tree lives.
        """
        actions: list[ModerationAction] = []

        # Example: if on cooldown, block message
        if self.is_on_cooldown(username):
            actions.append(
                ModerationAction(
                    action=ActionType.NONE,
                    reason="Message blocked: user on cooldown."
                )
            )
            return actions

        # Plug in rule modules here
        actions.extend(self._rule_bullshit_usage(username, message, chat_history))
        # actions.extend(self._rule_toxicity(username, message, chat_history))
        # actions.extend(self._rule_spam(username, message, chat_history))

        return actions