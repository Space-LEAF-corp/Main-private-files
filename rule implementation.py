import re

class EthicsEngine(EthicsEngine):  # extend previous
    def _mentions(self, message: str) -> list[str]:
        return re.findall(r'@([\w\-]+)', message)

    def _count_uses_against(self, username: str, target: str, chat_history: list) -> int:
        count = 0
        for sender, msg, _ in chat_history[-50:]:
            if sender == username and self.config.BAD_WORD in msg.lower() and f"@{target}" in msg:
                count += 1
        return count

    def _is_inappropriate_context(self, message: str) -> bool:
        lower_msg = message.lower()
        if self.config.BAD_WORD in lower_msg:
            return any(kw in lower_msg for kw in self.config.INAPPROPRIATE_KEYWORDS)
        return False

    def _rule_bullshit_usage(self, username: str, message: str, chat_history: list) -> list[ModerationAction]:
        actions = []
        lower_msg = message.lower()

        if self.config.BAD_WORD not in lower_msg:
            return actions

        # 1. Racial / inappropriate context → immediate ban (non‐weaponization)
        if self._is_inappropriate_context(message):
            actions.append(
                self.record_ban(
                    username,
                    self.config.BAN_DURATION_MINOR
                )
            )
            actions[-1].reason = "Inappropriate racial use of 'bullshit'."
            return actions

        # 2. Overuse against specific users → warnings → cooldown → ban
        mentions = self._mentions(message)
        for target in mentions:
            uses = self._count_uses_against(username, target, chat_history)

            if uses == self.config.WARNING_THRESHOLD:
                actions.append(
                    ModerationAction(
                        action=ActionType.WARN,
                        reason=f"Warning: overusing 'bullshit' toward @{target}."
                    )
                )
            elif uses == self.config.OVERUSE_THRESHOLD:
                # First step: cooldown, not ban
                actions.append(
                    self.record_cooldown(username, self.config.COOLDOWN_DURATION)
                )
                actions[-1].reason = f"Cooldown: repeated 'bullshit' toward @{target}."
            elif uses > self.config.OVERUSE_THRESHOLD:
                actions.append(
                    self.record_ban(
                        username,
                        self.config.BAN_DURATION_MAJOR
                    )
                )
                actions[-1].reason = f"Ban: persistent misuse of 'bullshit' toward @{target}."

        # 3. Appropriate flagging: "bullshit @user" → flag target for review
        if lower_msg.startswith(self.config.BAD_WORD) and mentions:
            target = mentions[0]
            actions.append(
                ModerationAction(
                    action=ActionType.FLAG,
                    reason=f"Flagged @{target} for review.",
                    target_user=target
                )
            )

        return actions